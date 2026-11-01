"""Offline GRCh37 -> GRCh38 point liftOver from the UCSC chain file.

Why local
---------
The Ensembl REST mapper answers one SNP per request at ~4 s each, i.e. hours
for 2 198 instruments, and kept timing out. The UCSC chain file (227 KB) is
parsed once; mapping is then instant and reproducible.

Chain convention (verified against the header sizes)
---------------------------------------------------
header: chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
In hg19/liftOver/hg19ToHg38.over.chain the **target (t) side is hg19**
(tSize 249250621 = hg19 chr1) and the **query (q) side is hg38**
(qSize 248956422), i.e. t is the source we lift *from*. Block lines
"size dt dq" are 0-based half-open. qStrand is '-' for inverted segments, in
which case qStart/qEnd are in reverse-complement space and a plus-strand
output needs qSize-1-q. Chain files are ordered by score, so when several
chains cover a base the earliest one wins.
"""
import bisect
import gzip
import os

CHAIN = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "data", "hg19ToHg38.over.chain.gz")


def _keep(name):
    return "_" not in name and name not in ("chrM", "chrUn")


class ChainLift:
    def __init__(self, path=CHAIN):
        # per chrom: (src_start0, src_end0, dst_start0, dst_strand, dst_size, rank)
        self.blocks = {}
        self._parse(path)
        self._keys = {}
        for c, v in self.blocks.items():
            v.sort(key=lambda b: b[0])
            self._keys[c] = [b[0] for b in v]
        self.n_blocks = sum(len(v) for v in self.blocks.values())

    def _parse(self, path):
        state = None
        rank = 0
        with gzip.open(path, "rt", errors="ignore") as fh:
            for line in fh:
                if line.startswith("chain"):
                    f = line.split()
                    t_name = f[2]
                    chrom = t_name if t_name.startswith("chr") else "chr" + t_name
                    rank += 1
                    if not _keep(chrom) or f[4] != "+":
                        state = None
                        continue
                    state = dict(chrom=chrom, src=int(f[5]), dst=int(f[10]),
                                 dst_strand=f[9], dst_size=int(f[8]), rank=rank)
                    continue
                if state is None:
                    continue
                f = line.split()
                if not f:
                    continue
                size = int(f[0])
                st = state
                self.blocks.setdefault(st["chrom"], []).append(
                    (st["src"], st["src"] + size, st["dst"],
                     st["dst_strand"], st["dst_size"], st["rank"]))
                if len(f) >= 3:
                    st["src"] += size + int(f[1])
                    st["dst"] += size + int(f[2])
                else:
                    state = None

    def map(self, chrom, pos1):
        """1-based GRCh37 position -> 1-based GRCh38 position, or None."""
        c = str(chrom)
        if not c.startswith("chr"):
            c = "chr" + c
        bl = self.blocks.get(c)
        if not bl:
            return None
        p0 = int(pos1) - 1
        i = bisect.bisect_right(self._keys[c], p0) - 1
        best = None
        j = i
        while j >= 0 and j > i - 300:
            s0, s1, d0, strand, dsize, rank = bl[j]
            if s1 > p0 >= s0:
                off = p0 - s0
                if strand == "+":
                    d = d0 + off
                else:
                    d = dsize - 1 - (d0 + off)
                if best is None or rank < best[0]:
                    best = (rank, d + 1)
            j -= 1
        return None if best is None else best[1]


if __name__ == "__main__":
    L = ChainLift()
    print("chroms: %d  blocks: %d" % (len(L.blocks), L.n_blocks))
    for c, p, expect in [("19", 726482, 726482), ("1", 113457000, 112914378),
                         ("7", 17276157, 17236533)]:
        got = L.map(c, p)
        print("%s:%-11d -> %-11s  ensembl=%s  %s"
              % (c, p, got, expect, "OK" if got == expect else "CHECK"))
    for c, p in [("19", 881161), ("6", 32600000), ("22", 22000000), ("17", 17600000)]:
        print("%s:%-11d -> %s" % (c, p, L.map(c, p)))
