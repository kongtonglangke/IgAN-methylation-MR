"""Minimal tabix/BGZF random-access reader over plain HTTPS Range requests.

Why this exists
---------------
EBI throttles bulk downloads from this host to ~13 KB/s per connection, so
pulling a 540-650 MB GWAS summary-statistics file takes hours. Almost every MR
lookup only needs a few hundred SNP positions, i.e. a few hundred 64 KB BGZF
blocks. This module fetches only those blocks.

Requirements: a remote `<file>.gz` plus its `<file>.gz.tbi` (tabix index).
EBI's GWAS Catalog "harmonised" files ship both.

Public API
----------
reader = TabixHTTP(url_gz, url_tbi, cache_dir="data/tbi_cache")
rows = reader.fetch([("1", 113286780), ("19", 726482), ...], chrom_col=0, pos_col=1)
    -> {("1", 113286780): [ [fields...], ... ]}

Notes
-----
* Coordinates are 1-based inclusive, matching both tabix and GWAS-SSF.
* Blocks are cached on disk so repeated runs cost nothing.
* Threaded: independent queries run in parallel; the block cache dedupes work.
"""
import gzip
import io
import os
import re
import ssl
import struct
import threading
import time
import urllib.request
import zlib
from concurrent.futures import ThreadPoolExecutor

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

MIN_SHIFT = 14
DEPTH = 5


def _range(url, start, end, timeout=120, retries=4):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=%d-%d" % (start, end)})
            with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
                return r.read()
        except Exception as ex:          # flaky FTP-over-HTTPS: retry with backoff
            last = ex
            time.sleep(2 + 3 * i)
    raise last


def reg2bins(beg, end):
    """htslib reg2bins for a 1-based inclusive region. Returns a set of bin ids."""
    if beg > end:
        return set()
    beg -= 1                      # to 0-based
    end -= 1
    s = MIN_SHIFT + DEPTH * 3
    bins = {0}
    t = 0
    for lvl in range(1, DEPTH + 1):
        s -= 3
        t += 1 << ((lvl - 1) * 3)
        b = t + (beg >> s)
        e = t + (end >> s)
        for b_ in range(b, e + 1):
            bins.add(b_)
    return bins


class TabixHTTP:
    def __init__(self, url_gz, url_tbi, cache_dir="data/tbi_cache", verbose=True):
        self.url = url_gz
        self.cache_dir = cache_dir
        self.verbose = verbose
        os.makedirs(cache_dir, exist_ok=True)
        self._lock = threading.Lock()
        self._blocks = {}
        self._idx_bytes = None
        self._load_index(url_tbi)

    # ---------- index ----------
    def _load_index(self, url_tbi):
        local = os.path.join(self.cache_dir, "index.tbi")
        if os.path.exists(local) and os.path.getsize(local) > 0:
            raw = open(local, "rb").read()
        else:
            raw = _range(url_tbi, 0, 1 << 31, timeout=300)
            open(local, "wb").write(raw)
        buf = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
        self._idx_bytes = buf
        self._parse_index(buf)

    def _parse_index(self, b):
        assert b[:4] == b"TBI\x01", "not a tabix index"
        off = 4
        (self.n_ref, self.fmt, self.col_seq, self.col_beg, self.col_end,
         self.meta, self.skip, l_nm) = struct.unpack_from("<8i", b, off)
        off += 32
        names = b[off:off + l_nm].split(b"\x00")
        off += l_nm
        self.names = [n.decode() for n in names if n]
        self.refs = {}
        for i in range(self.n_ref):
            (n_bin,) = struct.unpack_from("<i", b, off)
            off += 4
            bins = {}
            for _ in range(n_bin):
                (bin_id,) = struct.unpack_from("<I", b, off)
                (n_chunk,) = struct.unpack_from("<i", b, off + 4)
                off += 8
                ch = []
                for _ in range(n_chunk):
                    beg, end = struct.unpack_from("<QQ", b, off)
                    off += 16
                    ch.append((beg, end))
                bins[bin_id] = ch
            (n_intv,) = struct.unpack_from("<i", b, off)
            off += 4
            lin = []
            for _ in range(n_intv):
                (ioff,) = struct.unpack_from("<Q", b, off)
                off += 8
                lin.append(ioff)
            self.refs[self.names[i]] = {"bins": bins, "lin": lin}
        if self.verbose:
            tot = sum(len(c) for r in self.refs.values() for c in r["bins"].values())
            print("  index: %d refs, %d chunks" % (len(self.refs), tot), flush=True)

    # ---------- bgzf blocks ----------
    def _block(self, coffset):
        with self._lock:
            if coffset in self._blocks:
                return self._blocks[coffset]
        disk = os.path.join(self.cache_dir, "b_%d.bgz" % coffset)
        if os.path.exists(disk):
            blob = open(disk, "rb").read()
        else:
            head = _range(self.url, coffset, coffset + 1023)
            if head[:2] != b"\x1f\x8b":
                raise IOError("bad gzip magic at %d" % coffset)
            xlen = struct.unpack_from("<H", head, 10)[0]
            extra = head[12:12 + xlen]
            bsize = None
            i = 0
            while i + 4 <= len(extra):
                si = extra[i:i + 2]
                slen = struct.unpack_from("<H", extra, i + 2)[0]
                if si == b"BC":
                    bsize = struct.unpack_from("<H", extra, i + 4)[0]
                    break
                i += 4 + slen
            if bsize is None:
                raise IOError("no BC subfield at %d" % coffset)
            total = bsize + 1
            if len(head) >= total:
                blob = head[:total]
            else:
                blob = head + _range(self.url, coffset + len(head), coffset + total - 1)
                blob = blob[:total]
            open(disk, "wb").write(blob)
        dec = zlib.decompressobj(31).decompress(blob)
        out = (len(blob), dec)
        with self._lock:
            self._blocks[coffset] = out
        return out

    # ---------- query ----------
    def _query_lines(self, chrom, beg, end, max_blocks=60):
        """Return the tab-delimited lines overlapping [beg, end], reading only what
        is needed. The root bin (bin 0) spans the whole file, so the raw chunk list
        must be streamed and abandoned as soon as records pass `end`."""
        ref = self.refs.get(chrom)
        if ref is None:
            return []
        bins = reg2bins(beg, end)
        lin = ref["lin"]
        ki = beg >> MIN_SHIFT
        min_off = lin[ki] if lin and ki < len(lin) else 0
        chunks = []
        for b in bins:
            for cb, ce in ref["bins"].get(b, []):
                if ce > min_off:
                    chunks.append((max(cb, min_off), ce))
        if not chunks:
            return []
        chunks.sort()
        merged = [list(chunks[0])]
        for cb, ce in chunks[1:]:
            if cb <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], ce)
            else:
                merged.append([cb, ce])

        lines_out = []
        for cb, ce in merged:
            coff = cb >> 16
            uoff = cb & 0xFFFF
            end_coff = ce >> 16
            buf = b""
            first = True
            nblk = 0
            stop = False
            while True:
                nbytes, dec = self._block(coff)
                nblk += 1
                buf += dec[uoff:] if first else dec
                first = False
                parts = buf.split(b"\n")
                buf = parts.pop()
                for raw in parts:
                    if not raw or raw[:1] == b"#":
                        continue
                    f = raw.rstrip(b"\r").split(b"\t")
                    if len(f) < 2:
                        continue
                    try:
                        p = int(f[1])
                    except ValueError:
                        continue
                    if p > end:
                        stop = True
                        break
                    if p >= beg:
                        lines_out.append(raw.decode("utf-8", "ignore"))
                if stop:
                    break
                if coff + nbytes > end_coff or nblk >= max_blocks:
                    break
                coff += nbytes
                uoff = 0
        return lines_out

    def query(self, chrom, beg, end, max_blocks=60):
        """Public region query: returns list of tab-delimited field lists."""
        out = []
        for line in self._query_lines(str(chrom), int(beg), int(end), max_blocks=max_blocks):
            f = line.rstrip("\r").split("\t")
            if len(f) >= 2:
                out.append(f)
        return out

    def fetch(self, positions, nthreads=8, callback=None):
        """positions: iterable of (chrom, pos_1based). Returns dict key->list of field-lists."""
        want = {}
        for c, p in positions:
            want.setdefault((str(c), int(p)), []).append(True)

        def one(key):
            chrom, pos = key
            try:
                lines = self._query_lines(chrom, pos, pos)
            except Exception as ex:
                return (key, "ERR:" + type(ex).__name__ + ":" + str(ex)[:40])
            got = []
            for line in lines:
                f = line.rstrip("\r").split("\t")
                if len(f) < 2:
                    continue
                try:
                    if int(f[1]) != pos:
                        continue
                except ValueError:
                    continue
                c2 = f[0]
                if c2.lower().startswith("chr"):
                    c2 = c2[3:]
                if c2 != chrom:
                    continue
                got.append(f)
            return (key, got)

        res = {}
        done = 0
        with ThreadPoolExecutor(max_workers=nthreads) as ex:
            for key, val in ex.map(one, list(want.keys())):
                res[key] = val
                done += 1
                if callback and done % 50 == 0:
                    callback(done, len(want))
        return res
