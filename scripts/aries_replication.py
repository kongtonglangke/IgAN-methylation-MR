#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
aries_replication.py —— ARIES(mQTLdb) vs GoDMC 的独立队列方向一致性检验

为什么只能做"方向一致性"而不能做 coloc
------------------------------------
mQTLdb 的 MatrixEQTL 输出**只保留 FDR<0.05 的 CpG-SNP 配对**，没有未显著 SNP、
没有 varbeta 的完整区域信息 → coloc.abf 的前提（区域内逐 SNP 的 beta/SE）不成立。
因此本脚本只回答两个问题：
  Q1 在两边都检出的同一 SNP-CpG 配对上，效应方向是否一致？效应量是否相关？
  Q2 每个 CpG 的 GoDMC 索引 SNP 在 ARIES 中是否仍是同向显著 mQTL？（功效受限，阴性不可解读为"无效应"）

等位基因对齐
-----------
ARIES: A1 是效应等位（PLINK --freq 定义），beta 为秩正态化后甲基化的回归系数（SD 单位）
GoDMC: a1/a2，beta_a1 为 a1 的效应
两者同为 SD 单位 → 可直接比大小与方向。
互补链用 A<->T、C<->G 转换；A/T、C/G 这类无法用等位本身定链的，用 MAF 判向，仍不确定则丢弃。
"""
import os, io, json, math, glob, itertools, csv

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
COLOC = os.path.join(ROOT, "data", "coloc")

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
AMBIG = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

REGION = {
    "cg07209298": "CFH_CFHR_1q32", "cg11348583": "CFH_CFHR_1q32",
    "cg13441112": "CFH_CFHR_1q32", "cg13682187": "CFH_CFHR_1q32",
    "cg00705004": "DEFA1_8p23", "cg05969298": "DEFA1_8p23",
    "cg07297068": "DEFA1_8p23", "cg19304352": "DEFA1_8p23",
    "cg24860332": "DEFA1_8p23", "cg02787689": "DEFA1_8p23", "cg05787038": "DEFA1_8p23",
    "cg00539347": "HORMAD2_22q12", "cg27665648": "HORMAD2_22q12",
    "cg02846316": "ITGAM_ITGAX_16p11", "cg03100639": "ITGAM_ITGAX_16p11",
    "cg10451425": "ITGAM_ITGAX_16p11", "cg04018474": "ITGAM_ITGAX_16p11",
    "cg27535410": "MADCAM1_19p13", "cg23386659": "MADCAM1_19p13",
    "cg25599673": "ST6GAL1_3q27", "cg17448192": "ST6GAL1_3q27",
    "cg02671941": "TNFRSF13B_17p11",
}


def load_aries():
    """读 results/aries_mqtl_<cpg>.tsv，返回 dict cpg -> list of dict"""
    out = {}
    for p in sorted(glob.glob(os.path.join(RES, "aries_mqtl_cg*.tsv"))):
        cpg = os.path.basename(p)[len("aries_mqtl_"):-4]
        rows = []
        with io.open(p, encoding="utf-8") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for r in rd:
                try:
                    rows.append({
                        "tp": r["timepoint"], "snp": r["snp"],
                        "a1": r["a1"].upper(), "a2": r["a2"].upper(),
                        "maf": float(r["maf"]) if r["maf"] not in ("", "NA") else None,
                        "beta": float(r["beta"]), "t": float(r["tstat"]),
                        "p": float(r["pval"]),
                    })
                except (ValueError, KeyError):
                    continue
        if rows:
            out[cpg] = rows
    return out


def orient(ga1, ga2, gbeta, a1, a2, abeta, amaf, gfreq):
    """返回 (orientation, aligned_beta_aries)；orientation ∈ {same, flipped, None}"""
    if (a1, a2) == (ga1, ga2):
        return "same", abeta
    if (a1, a2) == (ga2, ga1):
        return "flipped", -abeta
    ca1, ca2 = COMP.get(a1), COMP.get(a2)
    if (ca1, ca2) == (ga1, ga2):
        return "same", abeta
    if (ca1, ca2) == (ga2, ga1):
        return "flipped", -abeta
    # 等位本身无法定链（A/T、C/G）→ 用频率判向
    if (a1, a2) in AMBIG and amaf is not None and gfreq is not None:
        if abs(amaf - gfreq) <= abs(amaf - (1 - gfreq)):
            return "same_freq", abeta
        return "flipped_freq", -abeta
    return None, None


def rank(v):
    """平均秩（含并列）"""
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def pearson(x, y):
    n = len(x)
    if n < 3:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def spearman(x, y):
    return pearson(rank(x), rank(y))


def binom_two_sided(k, n, p=0.5):
    """精确二项双侧 p；n 很大时用正态近似（含连续性校正）"""
    if n == 0:
        return float("nan")
    if n > 500:                       # 避免 comb 溢出，且精确解在此规模已无必要
        mu, sd = n * p, math.sqrt(n * p * (1 - p))
        z = (abs(k - mu) - 0.5) / sd
        return max(0.0, min(1.0, math.erfc(z / math.sqrt(2))))

    def logpmf(i):
        return (math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                + i * math.log(p) + (n - i) * math.log(1 - p))
    obs = logpmf(k)
    tot = 0.0
    for i in range(n + 1):
        v = logpmf(i)
        if v <= obs + 1e-9:
            tot += math.exp(v)
    return min(1.0, tot)


def main():
    aries = load_aries()
    pairs, summ = [], []
    for cpg, arows in sorted(aries.items()):
        gp = os.path.join(COLOC, "mqtl_%s.json" % cpg)
        if not os.path.exists(gp):
            continue
        godmc = {x["rsid"]: x for x in json.load(open(gp)) if x.get("rsid")}
        best = min(godmc.values(), key=lambda x: x["pval_mre"]) if godmc else None
        idx_rs = best["rsid"] if best else None
        seen, local = set(), []
        for a in arows:
            g = godmc.get(a["snp"])
            if not g:
                continue
            o, ab = orient(g["a1"].upper(), g["a2"].upper(), g["beta_a1"],
                           a["a1"], a["a2"], a["beta"], a["maf"], g.get("freq_a1"))
            if o is None:
                continue
            key = (a["tp"], a["snp"])
            if key in seen:
                continue
            seen.add(key)
            conc = (ab > 0) == (g["beta_a1"] > 0)
            local.append({"cpg": cpg, "region": REGION.get(cpg, "NA"),
                          "timepoint": a["tp"], "snp": a["snp"],
                          "aries_a1": a["a1"], "aries_a2": a["a2"], "aries_maf": a["maf"],
                          "godmc_a1": g["a1"], "godmc_a2": g["a2"], "godmc_freq": g.get("freq_a1"),
                          "orient": o, "aries_beta_aligned": ab, "aries_t": a["t"], "aries_p": a["p"],
                          "godmc_beta": g["beta_a1"], "godmc_p_mre": g["pval_mre"],
                          "godmc_I2": g.get("hetisq"), "n_studies": g.get("num_studies"),
                          "concordant": int(conc), "is_index_snp": int(a["snp"] == idx_rs)})
        if not local:
            continue
        pairs.extend(local)
        n = len(local)
        k = sum(x["concordant"] for x in local)
        ab = [x["aries_beta_aligned"] for x in local]
        gb = [x["godmc_beta"] for x in local]
        idxrows = [x for x in local if x["is_index_snp"]]
        summ.append({
            "cpg": cpg, "region": REGION.get(cpg, "NA"), "n_pairs": n, "n_conc": k,
            "pct_conc": 100.0 * k / n, "binom_p": binom_two_sided(k, n),
            "spearman": spearman(ab, gb), "pearson": pearson(ab, gb),
            "index_snp": idx_rs or "",
            "index_in_aries": int(bool(idxrows)),
            "index_aries_beta": idxrows[0]["aries_beta_aligned"] if idxrows else "",
            "index_aries_p": idxrows[0]["aries_p"] if idxrows else "",
            "index_godmc_beta": idxrows[0]["godmc_beta"] if idxrows else "",
            "index_same_dir": (int((idxrows[0]["aries_beta_aligned"] > 0) == (idxrows[0]["godmc_beta"] > 0))
                               if idxrows else ""),
        })

    # 写 pair 表
    pp = os.path.join(RES, "aries_vs_godmc_pairs.tsv")
    cols = ["cpg", "region", "timepoint", "snp", "aries_a1", "aries_a2", "aries_maf",
            "godmc_a1", "godmc_a2", "godmc_freq", "orient", "aries_beta_aligned", "aries_t",
            "aries_p", "godmc_beta", "godmc_p_mre", "godmc_I2", "n_studies",
            "concordant", "is_index_snp"]
    with io.open(pp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in pairs:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    sp = os.path.join(RES, "aries_replication_summary.tsv")
    scols = ["cpg", "region", "n_pairs", "n_conc", "pct_conc", "binom_p", "spearman",
             "pearson", "index_snp", "index_in_aries", "index_aries_beta", "index_aries_p",
             "index_godmc_beta", "index_same_dir"]
    with io.open(sp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(scols) + "\n")
        for r in summ:
            fh.write("\t".join(("%.4g" % v if isinstance(v, float) else str(v)) for v in
                               [r[c] for c in scols]) + "\n")

    # 总体
    tot_n = len(pairs)
    tot_k = sum(x["concordant"] for x in pairs)
    print("CpG 数: %d  配对总数: %d  同向: %d (%.1f%%)  二项检验 p=%.3g"
          % (len(summ), tot_n, tot_k, 100.0 * tot_k / max(tot_n, 1),
             binom_two_sided(tot_k, tot_n)))
    print("整体 Spearman(beta): %.3f" % spearman([x["aries_beta_aligned"] for x in pairs],
                                                 [x["godmc_beta"] for x in pairs]))
    print("\n%-12s %-18s %6s %6s %7s %8s %8s" % ("cpg", "region", "n", "conc", "%", "rho", "binom_p"))
    for r in summ:
        print("%-12s %-18s %6d %6d %7.1f %8.3f %8.3g"
              % (r["cpg"], r["region"], r["n_pairs"], r["n_conc"], r["pct_conc"],
                 r["spearman"], r["binom_p"]))
    print("\n写出 %s / %s" % (pp, sp))


if __name__ == "__main__":
    main()
