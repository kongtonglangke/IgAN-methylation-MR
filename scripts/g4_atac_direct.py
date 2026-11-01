# -*- coding: utf-8 -*-
"""G4: 直读 ChIP-Atlas ATAC 的指定区域

tabix chunk 偏移为虚拟偏移 voffset=(cOffset<<16)|uOffset，取 cOffset 即 bgzf 块边界。
只读取覆盖目标区域的 bgzf 块，避免下载 59.5 GB 主文件。
"""
import gzip
import io
import os
import ssl
import struct
import sys
import time
import urllib.request
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tabix_http import reg2bins

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G4 = os.path.join(ROOT, "data", "g4")
RES = os.path.join(ROOT, "results")

GZ = "https://chip-atlas.dbcls.jp/data/hg38/assembled/ATC.ALL.05.AllAg.AllCell.bed.gz"
FILE_SIZE = 59455721996  # 来自 HEAD Content-Length


def load_index(path):
    raw = open(path, "rb").read()
    b = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
    assert b[:4] == b"TBI\x01"
    off = 4
    (n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm) = struct.unpack_from("<8i", b, off)
    off += 32
    names = b[off:off + l_nm].split(b"\x00")
    off += l_nm
    names = [n.decode() for n in names if n]
    refs = {}
    for i in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", b, off)
        off += 4
        bins = {}
        for _ in range(n_bin):
            (bid,) = struct.unpack_from("<I", b, off)
            (nch,) = struct.unpack_from("<i", b, off + 4)
            off += 8
            ch = []
            for _ in range(nch):
                bg, en = struct.unpack_from("<QQ", b, off)
                off += 16
                ch.append((bg, en))
            bins[bid] = ch
        (n_intv,) = struct.unpack_from("<i", b, off)
        off += 4
        lin = []
        for _ in range(n_intv):
            (ioff,) = struct.unpack_from("<Q", b, off)
            off += 8
            lin.append(ioff)
        refs[names[i]] = {"bins": bins, "lin": lin}
    return refs


def plan_chunks(ref, beg, end):
    """返回需读取的字节区间列表 [(bstart,bend)] 与总字节数"""
    b0 = beg - 1
    bins = reg2bins(b0, end)
    lin = ref["lin"]
    li = b0 >> 14
    min_voff = lin[li] if li < len(lin) else (lin[-1] if lin else 0)
    segs = []
    for x in bins:
        for cs, ce in ref["bins"].get(x, []):
            if ce > min_voff:
                segs.append((cs >> 16, ce >> 16))
    if not segs:
        return [], 0
    segs.sort()
    merged = [list(segs[0])]
    for s, e in segs[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(m) for m in merged], sum(e - s for s, e in merged)


def fetch(url, start, end, retries=4):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=%d-%d" % (start, end - 1)})
            with urllib.request.urlopen(req, timeout=180, context=CTX) as r:
                return r.read()
        except Exception as ex:
            last = ex
            time.sleep(3 * (i + 1))
    raise last


def gunzip_all(data):
    out = []
    rest = data
    while rest:
        d = zlib.decompressobj(31)
        try:
            out.append(d.decompress(rest))
        except zlib.error:
            break
        rest = d.unused_data
        if not rest:
            break
    return b"".join(out)


def query_region(refs, chrom, beg, end, verbose=True):
    ref = refs[chrom]
    segs, total = plan_chunks(ref, beg, end)
    if verbose:
        print("  chunks: %d  bytes: %d (%.2f MB)" % (len(segs), total, total / 1048576.0), flush=True)
    hits = []
    for k, (s, e) in enumerate(segs):
        e = min(e, FILE_SIZE)
        if s >= FILE_SIZE or e <= s:
            if verbose:
                print("    seg %d/%d skipped (out of file bounds)" % (k + 1, len(segs)), flush=True)
            continue
        t = time.time()
        raw = fetch(GZ, s, e)
        txt = gunzip_all(raw)
        for line in txt.decode("utf-8", "replace").splitlines():
            f = line.split("\t")
            if len(f) < 3:
                continue
            if f[0] != chrom:
                continue
            try:
                ps, pe = int(f[1]), int(f[2])
            except ValueError:
                continue
            if pe > beg and ps < end:
                hits.append(f)
        if verbose:
            print("    seg %d/%d %d bytes -> %d hits  %.1fs" %
                  (k + 1, len(segs), len(raw), len(hits), time.time() - t), flush=True)
    return hits


REGIONS = [
    ("DEFA1_8p23", "chr8", 6965292, 6975292, 6970292, "cg00705004"),
    ("HORMAD2_22q12", "chr22", 30191307, 30201307, 30196307, "cg00539347"),
]


def main():
    refs = load_index(os.path.join(G4, "tabix_cache", "index.tbi"))
    print("index: %d refs" % len(refs), flush=True)
    rows = []
    for label, chrom, s, e, pos, cpg in REGIONS:
        print("\n== %s %s:%d-%d  (CpG %d)" % (label, chrom, s, e, pos), flush=True)
        hits = query_region(refs, chrom, s, e)
        print("  peaks in window: %d" % len(hits), flush=True)
        if hits:
            print("  example: %s" % ("\t".join(hits[0][:6])), flush=True)
        for f in hits:
            ps, pe = int(f[1]), int(f[2])
            inside = ps <= pos < pe
            dist = 0 if inside else min(abs(pos - ps), abs(pos - pe))
            rows.append((label, cpg, chrom, ps, pe, pe - ps,
                         f[3] if len(f) > 3 else "", f[4] if len(f) > 4 else "",
                         "YES" if inside else "no", dist))
        ins = [r for r in rows if r[0] == label and r[8] == "YES"]
        print("  CpG inside peak: %s" % ("YES (%d peaks)" % len(ins) if ins else "no"), flush=True)
        if not ins and hits:
            near = min([r for r in rows if r[0] == label], key=lambda r: r[9])
            print("  nearest peak edge: %d bp" % near[9], flush=True)

    outp = os.path.join(RES, "g4_atac_peaks.tsv")
    with io.open(outp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(["region", "cpg", "chrom", "start", "end", "length_bp",
                            "name", "score", "cpg_inside", "dist_to_cpg_bp"]) + "\n")
        for r in rows:
            fh.write("\t".join(str(v) for v in r) + "\n")
    print("\n写出 %s (%d 行)" % (outp, len(rows)))


if __name__ == "__main__":
    main()
