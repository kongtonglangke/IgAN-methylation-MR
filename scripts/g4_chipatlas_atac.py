# -*- coding: utf-8 -*-
"""G4: ChIP-Atlas hg38 ATAC-seq (ATC.ALL.05.AllAg.AllCell) 区域随机读取

使用已下载的 tabix 索引 + HTTP Range，仅取目标区域，避免下载 59.5 GB 主文件。
tabix 坐标为 1-based inclusive。
"""
import io
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tabix_http import TabixHTTP

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G4 = os.path.join(ROOT, "data", "g4")
RES = os.path.join(ROOT, "results")
CACHE = os.path.join(G4, "tabix_cache")

GZ = "https://chip-atlas.dbcls.jp/data/hg38/assembled/ATC.ALL.05.AllAg.AllCell.bed.gz"

REGIONS = [
    ("DEFA1_8p23", "chr8", 6960000, 6982000, 6970292, "cg00705004"),
    ("HORMAD2_22q12", "chr22", 30188000, 30205000, 30196307, "cg00539347"),
]


def main():
    os.makedirs(CACHE, exist_ok=True)
    local_idx = os.path.join(CACHE, "index.tbi")
    src = os.path.join(G4, "atac_hg38.tbi")
    if not (os.path.exists(local_idx) and os.path.getsize(local_idx) > 1000000):
        shutil.copyfile(src, local_idx)
        print("copy index ->", local_idx, os.path.getsize(local_idx), flush=True)

    reader = TabixHTTP(GZ, "https://chip-atlas.dbcls.jp/data/hg38/assembled/ATC.ALL.05.AllAg.AllCell.bed.gz.tbi",
                       cache_dir=CACHE, verbose=True)

    rows = []
    for label, chrom, s, e, pos, cpg in REGIONS:
        print("\n== %s %s:%d-%d (CpG %d)" % (label, chrom, s, e, pos), flush=True)
        recs = reader.query(chrom, s + 1, e, max_blocks=200)
        print("  peaks returned: %d" % len(recs), flush=True)
        if recs:
            print("  first record fields: %s" % (recs[0][:8],), flush=True)
        for f in recs:
            try:
                ps, pe = int(f[1]), int(f[2])
            except (ValueError, IndexError):
                continue
            inside = ps <= pos < pe
            dist = 0 if inside else min(abs(pos - ps), abs(pos - pe))
            rows.append((label, cpg, chrom, ps, pe, pe - ps,
                         f[3] if len(f) > 3 else "", f[4] if len(f) > 4 else "",
                         "YES" if inside else "no", dist))
        ins = [r for r in rows if r[0] == label and r[8] == "YES"]
        if ins:
            print("  CpG INSIDE %d ATAC peak(s)" % len(ins), flush=True)
        elif rows:
            sub = [r for r in rows if r[0] == label]
            if sub:
                near = min(sub, key=lambda r: r[9])
                print("  CpG not in peak; nearest edge %d bp" % near[9], flush=True)

    outp = os.path.join(RES, "g4_atac_peaks.tsv")
    with io.open(outp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(["region", "cpg", "chrom", "start", "end", "length_bp",
                            "name", "score", "cpg_inside", "dist_to_cpg_bp"]) + "\n")
        for r in rows:
            fh.write("\t".join(str(v) for v in r) + "\n")
    print("\n写出 %s (%d 行)" % (outp, len(rows)))


if __name__ == "__main__":
    main()
