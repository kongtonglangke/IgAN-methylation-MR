# -*- coding: utf-8 -*-
"""G4: ENCODE cCRE 注释 + 区域调控密度对比

输入: data/g4/ccre_chr8_wide.json (chr8 110kb), ccre_chr22.json (chr22 40kb)
输出: results/g4_ccre_annotation.tsv, results/g4_ccre_density.tsv
"""
import io
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G4 = os.path.join(ROOT, "data", "g4")
RES = os.path.join(ROOT, "results")

# UCSC encodeCcreCombined 的 reserved 字段为 RGB；映射经 TSS 位置验证：
# 唯一 PLS 正对 LOC105372988 TSS(30,207,194)；pELS 聚集在 TSS 2kb 邻域
RGB2TYPE = {
    "255,0,0": "PLS",              # promoter-like
    "255,167,0": "pELS",           # proximal enhancer-like
    "255,205,0": "dELS",           # distal enhancer-like
    "0,176,240": "CTCF-only",
    "0,176,80": "DNase-H3K4me3",
}

TARGETS = [
    # (label, file, chrom, window_start, window_end, cpg_pos, cpg_id, gene_note)
    ("DEFA1_8p23", "ccre_chr8_wide.json", "chr8", 6920000, 7030000, 6970292,
     "cg00705004", "DEFA10P TSS-1179 / DEFA1 TSS+9800"),
    ("HORMAD2_22q12", "ccre_chr22.json", "chr22", 30176000, 30216000, 30196307,
     "cg00539347", "LOC105372988 内含子 / HORMAD2 TSS-115844"),
]

# 基因 TSS（GRCh38, ncbiRefSeqCurated）
GENE_TSS = {
    "chr8": [("DEFA10P", 6969113, "-"), ("DEFA1", 6980092, "-"),
             ("DEFT1P", 6989721, "-"), ("DEFA1B", 6999198, "-")],
    "chr22": [("LOC105372988", 30207194, "-"), ("HORMAD2", 30080463, "+")],
}


def classify(x):
    return RGB2TYPE.get(x["reserved"], x["reserved"])


def main():
    rows = []
    dens = []
    for label, fn, chrom, ws, we, pos, cpg, note in TARGETS:
        path = os.path.join(G4, fn)
        d = json.load(io.open(path, encoding="utf-8"))["encodeCcreCombined"]
        span = we - ws
        n = len(d)
        dens.append((label, chrom, ws, we, span, n,
                     round(n / (span / 1000.0), 3), cpg, note))
        for x in d:
            mid = (x["chromStart"] + x["chromEnd"]) // 2
            dist = mid - pos
            inside = x["chromStart"] <= pos < x["chromEnd"]
            nearest_edge = 0
            if not inside:
                nearest_edge = min(abs(pos - x["chromStart"]), abs(pos - x["chromEnd"]))
            # 到最近基因 TSS 的距离
            tss_best, tss_gene, tss_d = None, "", None
            for g, t, st in GENE_TSS.get(chrom, []):
                dd = abs(t - mid)
                if tss_d is None or dd < tss_d:
                    tss_d, tss_gene, tss_best = dd, g, t
            rows.append((label, cpg, chrom, x["name"], x["chromStart"], x["chromEnd"],
                         classify(x), x.get("score", ""), x["chromEnd"] - x["chromStart"],
                         dist, "YES" if inside else "no", nearest_edge,
                         tss_gene, tss_best, tss_d))

    outp = os.path.join(RES, "g4_ccre_annotation.tsv")
    with io.open(outp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(["region", "cpg", "chrom", "ccre_id", "start", "end",
                            "ccre_class", "encode_score", "length_bp",
                            "dist_to_cpg_bp", "cpg_inside", "dist_to_nearest_edge_bp",
                            "nearest_gene", "nearest_tss", "dist_to_tss_bp"]) + "\n")
        for r in rows:
            fh.write("\t".join(str(v) for v in r) + "\n")

    outd = os.path.join(RES, "g4_ccre_density.tsv")
    with io.open(outd, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(["region", "chrom", "window_start", "window_end",
                            "window_bp", "n_ccre", "ccre_per_kb", "cpg", "gene_note"]) + "\n")
        for r in dens:
            fh.write("\t".join(str(v) for v in r) + "\n")

    # ---- 控制台摘要 ----
    print("写出 %s (%d 行)" % (outp, len(rows)))
    print("写出 %s (%d 行)" % (outd, len(dens)))
    print()
    # 全基因组背景密度
    print("全基因组背景: ~2.35M cCRE / 3.09 Gb = 0.76 cCRE/kb (ENCODE Registry V4)")
    for label, chrom, ws, we, span, n, perkb, cpg, note in dens:
        print("  %-14s %s:%d-%d  %d kb  n=%d  %.3f/kb" % (label, chrom, ws, we, span // 1000, n, perkb))
    print()

    for label, fn, chrom, ws, we, pos, cpg, note in TARGETS:
        sub = [r for r in rows if r[0] == label]
        inside = [r for r in sub if r[10] == "YES"]
        print("== %s / %s (%s:%d)" % (label, cpg, chrom, pos))
        if inside:
            for r in inside:
                print("   INSIDE cCRE: %s %s %d-%d score=%s" % (r[3], r[6], r[4], r[5], r[7]))
        else:
            near = min(sub, key=lambda r: r[11])
            print("   NOT in any cCRE; nearest edge %d bp -> %s (%s) %d-%d" %
                  (near[11], near[3], near[6], near[4], near[5]))
        import collections
        print("   class counts:", dict(collections.Counter(r[6] for r in sub)))
        # 到最近基因 TSS 是否有 cCRE
        for g, t, st in GENE_TSS.get(chrom, []):
            cand = [r for r in sub if abs(((r[4] + r[5]) // 2) - t) <= 2000]
            print("   TSS of %-14s (%d): cCRE within 2kb = %d" % (g, t, len(cand)))
        print()


if __name__ == "__main__":
    main()
