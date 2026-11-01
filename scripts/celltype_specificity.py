#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
celltype_specificity.py —— 用 GSE35069（Reinius 2012, 纯化血细胞 450K, n=60）
判断 19 个 coloc CpG 的甲基化在哪一类血细胞里最高/最低。

用途
----
1. 若某 CpG 在粒细胞/中性粒细胞中显著高于其他细胞 → 该位点的甲基化信号更可能是
   髓系/中性粒细胞来源，而不是 B 细胞/浆细胞来源，"B 细胞表观遗传枢纽"的叙事要改写。
2. 为"细胞类型归属"提供可复算的实证，而不是靠文献推测。

输入: data/geo/GSE35069_series_matrix.txt.gz
输出: results/celltype_specificity.tsv
"""
import os, io, gzip, math, csv

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MTX = os.path.join(ROOT, "data", "geo", "GSE35069_series_matrix.txt.gz")
OUT = os.path.join(ROOT, "results", "celltype_specificity.tsv")

TARGET = ["cg07209298", "cg11348583", "cg13441112", "cg13682187",
          "cg00705004", "cg05969298", "cg07297068", "cg19304352",
          "cg24860332", "cg02787689", "cg05787038",
          "cg00539347", "cg27665648",
          "cg02846316", "cg03100639", "cg10451425", "cg04018474",
          "cg27535410", "cg23386659", "cg25599673", "cg17448192", "cg02671941"]

# GSE35069 实际样本标题为 WB / PBMC / Gran / CD4+ / CD8+ / CD14+ / CD19+ / CD56+ / Neu / Eos
# （各 6 例，共 60 例）；直接按前缀映射，别用模糊关键字
PREFIX = [
    ("CD19", "B_cell"), ("CD14", "Monocyte"), ("CD4", "CD4_T"), ("CD8", "CD8_T"),
    ("CD56", "NK"), ("Neu", "Neutrophil"), ("Gran", "Granulocyte"),
    ("Eos", "Eosinophil"), ("PBMC", "PBMC"), ("WB", "WholeBlood"),
]


def classify(title):
    t = title.split("_")[0].upper()
    for pre, grp in PREFIX:
        if t.startswith(pre.upper()):
            return grp
    return None


def main():
    if not os.path.exists(MTX):
        raise SystemExit("缺 %s" % MTX)
    titles, accs, data = [], [], {}
    with gzip.open(MTX, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accs = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                break
        rd = csv.reader(fh, delimiter="\t", quotechar='"')
        hdr = next(rd)
        idx = {g: i for i, g in enumerate(hdr)}
        want = set(TARGET)
        for row in rd:
            if not row or row[0].startswith("!series_matrix_table_end"):
                break
            cg = row[0].strip('"')
            if cg in want:
                vals = []
                for g in accs:
                    j = idx.get(g)
                    try:
                        vals.append(float(row[j]) if j is not None and row[j] not in ("", "null") else None)
                    except ValueError:
                        vals.append(None)
                data[cg] = vals

    groups = {}
    for i, t in enumerate(titles):
        g = classify(t)
        if g:
            groups.setdefault(g, []).append(i)
    order = ["Neutrophil", "Granulocyte", "Eosinophil", "Monocyte", "B_cell",
             "CD4_T", "CD8_T", "NK", "PBMC", "WholeBlood"]
    order = [g for g in order if g in groups]

    rows = []
    for cg in TARGET:
        if cg not in data:
            rows.append({"cpg": cg, "status": "NOT_ON_ARRAY"})
            continue
        v = data[cg]
        means = {}
        for g in order:
            xs = [v[i] for i in groups[g] if v[i] is not None]
            means[g] = sum(xs) / len(xs) if xs else float("nan")
        vals = [means[g] for g in order if means[g] == means[g]]
        if not vals:
            rows.append({"cpg": cg, "status": "NO_DATA"})
            continue
        hi = max(order, key=lambda g: means[g] if means[g] == means[g] else -9)
        lo = min(order, key=lambda g: means[g] if means[g] == means[g] else 9)
        rng = means[hi] - means[lo]
        my = sum(vals) / len(vals)
        sd = math.sqrt(sum((x - my) ** 2 for x in vals) / (len(vals) - 1)) if len(vals) > 1 else 0.0
        r = {"cpg": cg, "status": "OK", "max_cell": hi, "max_beta": means[hi],
             "min_cell": lo, "min_beta": means[lo], "range": rng,
             "sd_across_celltypes": sd, "cv": (sd / my if my else float("nan"))}
        for g in order:
            r[g] = means[g]
        rows.append(r)

    cols = ["cpg", "status"] + order + ["max_cell", "max_beta", "min_cell", "min_beta",
                                       "range", "sd_across_celltypes", "cv"]
    with io.open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(("%.4f" % r[c] if isinstance(r.get(c), float) else str(r.get(c, "")))
                               for c in cols) + "\n")

    print("样本分类:", {g: len(groups[g]) for g in order})
    print("%-12s %-12s %7s %7s %7s  %s" % ("cpg", "status", "max", "min", "range", "最高细胞"))
    for r in rows:
        if r["status"] != "OK":
            print("%-12s %-12s" % (r["cpg"], r["status"]))
            continue
        print("%-12s %-12s %7.3f %7.3f %7.3f  %s" % (r["cpg"], r["status"], r["max_beta"],
                                                     r["min_beta"], r["range"], r["max_cell"]))
    print("写出", OUT)


if __name__ == "__main__":
    main()
