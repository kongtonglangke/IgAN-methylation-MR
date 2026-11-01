#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
g1_feasibility.py —— 「菌群 → CpG」效应矩阵的可行性预跑

背景
----
Phase 1 的两步中介 MR 卡在"缺菌群→CpG 效应矩阵"。理论上可用两样本 MR 补：
  exposure = MiBioGen 菌群丰度（工具 SNP）
  outcome  = GoDMC 甲基化（cis-mQTL 区域）
但先要回答：**MiBioGen 的工具 SNP 在目标区域内到底有没有覆盖？**
若每个 taxa 在这些区域只有 0–2 个工具 SNP，整套分析直接判死，不必跑。

输入
----
_g0/MBG_allHits.txt  = MiBioGen MBG.allHits.p1e4.txt（全部 p<1e-4 关联）
    列: bac, chr, bp(GRCh37), rsID, ref.allele, eff.allele, beta, SE, Z, P, N, Ncohorts

输出
----
results/g1_microbiota_region_coverage.tsv   每个区域×taxa 的工具 SNP 数
results/g1_microbiota_region_summary.tsv    每个区域的汇总
"""
import io, os, csv, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = "F:/longxiawenjian/qiguanzhou/_g0/MBG_allHits.txt"
RES = os.path.join(ROOT, "results")

# GRCh37；窗口取 CpG 位置 ±500 kb（与 GoDMC cis 窗口一致）
REGIONS = [
    ("DEFA1_8p23",       8,  6_827_814),
    ("HORMAD2_22q12",   22, 30_143_693),
    ("CFH_CFHR_1q32",    1, 196_717_788),
    ("ITGAM_ITGAX_16p11", 16, 31_346_439),
    ("MADCAM1_19p13",   19,    843_692),
    ("ST6GAL1_3q27",     3, 187_020_633),
    ("TNFRSF13B_17p11", 17, 16_948_136),
]
WIN = 500_000
P_INST = 1e-5          # 工具变量阈值（MiBioGen 通用）


def main():
    if not os.path.exists(SRC):
        raise SystemExit("缺 %s" % SRC)

    # 逐区域收集
    hit = collections.defaultdict(list)      # (region, taxon) -> [(rsid, beta, se, p)]
    all_taxa = set()
    with io.open(SRC, encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        for row in rd:
            if len(row) < 12:
                continue
            taxon = row[0].strip('"')
            try:
                chrom = int(row[1]); pos = int(row[2])
                beta = float(row[6]); se = float(row[7]); p = float(row[9])
            except ValueError:
                continue
            all_taxa.add(taxon)
            if p > P_INST or se <= 0:
                continue
            for name, c, anchor in REGIONS:
                if chrom == c and abs(pos - anchor) <= WIN:
                    hit[(name, taxon)].append((row[3].strip('"'), beta, se, p))

    rows = []
    for (region, taxon), lst in sorted(hit.items()):
        f = [(b / s) ** 2 for _, b, s, _ in lst]
        rows.append({"region": region, "taxon": taxon, "n_snp": len(lst),
                     "min_p": min(x[3] for x in lst), "median_F": sorted(f)[len(f) // 2],
                     "max_F": max(f), "rsids": ",".join(x[0] for x in lst[:8])})
    outp = os.path.join(RES, "g1_microbiota_region_coverage.tsv")
    cols = ["region", "taxon", "n_snp", "min_p", "median_F", "max_F", "rsids"]
    with io.open(outp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(("%.4g" % r[c] if isinstance(r[c], float) else str(r[c]))
                               for c in cols) + "\n")

    print("MiBioGen taxa 总数: %d ；p<1e-5 的工具 SNP 命中区域:" % len(all_taxa))
    print("%-20s %8s %8s %10s %10s" % ("region", "n_taxa", "n_snp", "max_snp/taxa", "min_p"))
    srows = []
    for name, c, anchor in REGIONS:
        rel = [r for r in rows if r["region"] == name]
        n_snp = sum(r["n_snp"] for r in rel)
        mx = max([r["n_snp"] for r in rel], default=0)
        mp = min([r["min_p"] for r in rel], default=float("nan"))
        print("%-20s %8d %8d %10d %10.2e" % (name, len(rel), n_snp, mx, mp))
        srows.append({"region": name, "n_taxa": len(rel), "n_snp": n_snp,
                      "max_snp_per_taxon": mx, "min_p": mp})
    sp = os.path.join(RES, "g1_microbiota_region_summary.tsv")
    with io.open(sp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("region\tn_taxa\tn_snp\tmax_snp_per_taxon\tmin_p\n")
        for r in srows:
            fh.write("%s\t%d\t%d\t%d\t%.4g\n" % (r["region"], r["n_taxa"], r["n_snp"],
                                                 r["max_snp_per_taxon"], r["min_p"]))
    print("\n写出", outp, "与", sp)


if __name__ == "__main__":
    main()
