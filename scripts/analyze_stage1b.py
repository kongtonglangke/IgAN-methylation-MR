"""Summarise stage-1b (CpG methylation -> IgA nephropathy).

Reports, in order: instrument coverage in the outcome GWAS, the global
p-value distribution (calibration), locus-level and CpG-level hits, and the
detectable-effect floor implied by the observed standard errors.
"""
import csv
import math
from collections import defaultdict

import numpy as np


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = o[rank]
        prev = min(prev, p[i] * n / (rank + 1))
        q[i] = prev
    return q


def ks_uniform(p):
    p = np.sort(np.asarray(p, float))
    n = len(p)
    d = max(np.max(np.arange(1, n + 1) / n - p), np.max(p - np.arange(0, n) / n))
    lam = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * d
    pv = 2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam) for k in range(1, 101))
    return d, min(max(pv, 0.0), 1.0)


rows = [r for r in csv.DictReader(open("results/mr_stage1b_cpg_igan.tsv"), delimiter="\t")]
from collections import Counter
print("total CpGs screened : %d" % len(rows))
print("status              : %s" % dict(Counter(r["status"] for r in rows)))
ok = [r for r in rows if r["status"] == "ok"]
for r in ok:
    r["_p"] = float(r["pval"])
    r["_b"] = float(r["b"])
    r["_se"] = float(r["se"])
print("analysable          : %d CpGs / %d genes"
      % (len(ok), len({r["gene"] for r in ok})))

p = np.array([r["_p"] for r in ok])
n = len(p)
print("\n--- global calibration ---")
print("p<0.05 observed/expected : %d / %.1f" % ((p < 0.05).sum(), 0.05 * n))
print("p<0.01 observed/expected : %d / %.1f" % ((p < 0.01).sum(), 0.01 * n))
print("p<1e-3 observed/expected : %d / %.1f" % ((p < 1e-3).sum(), 1e-3 * n))
z = None
from scipy import stats
chi = stats.chi2.isf(p, 1)
lam = np.median(chi) / 0.4549
d, ksp = ks_uniform(p)
print("KS uniformity            : D=%.4f  p=%.3f" % (d, ksp))
print("lambda_GC                : %.3f" % lam)

q = bh(p)
nb = 0.05 / n
print("\n--- multiplicity ---")
print("Bonferroni threshold (n=%d): %.2e" % (n, nb))
print("CpGs passing Bonferroni    : %d" % (p < nb).sum())
print("CpGs passing FDR<0.05      : %d" % (q < 0.05).sum())
print("CpGs passing FDR<0.10      : %d" % (q < 0.10).sum())

print("\n--- top 15 CpGs by p ---")
print("%-12s %-10s %-12s %7s %9s %9s %10s" %
      ("cpg", "gene", "rsid", "F", "beta", "se", "p"))
for i in np.argsort(p)[:15]:
    r = ok[i]
    print("%-12s %-10s %-12s %7s %9.3f %9.3f %10.2e"
          % (r["cpg"], r["gene"], r["rsid"], r["lead_F"], r["_b"], r["_se"], r["_p"]))

print("\n--- locus-level (min-p per gene, Šidák within locus) ---")
byg = defaultdict(list)
for r in ok:
    byg[r["gene"]].append(r)
loci = []
for g, rs in byg.items():
    rs.sort(key=lambda x: x["_p"])
    k = len(rs)
    pmin = rs[0]["_p"]
    p_sid = 1 - (1 - pmin) ** k
    loci.append((g, k, pmin, p_sid, rs[0]["cpg"], rs[0]["_b"]))
loci.sort(key=lambda x: x[2])
print("%-10s %5s %11s %11s %-12s %8s" % ("gene", "nCpG", "p_min", "p_sidak", "cpg", "beta"))
for g, k, pmin, ps, cpg, b in loci:
    print("%-10s %5d %11.2e %11.2e %-12s %8.3f" % (g, k, pmin, ps, cpg, b))

print("\n--- detectable effect floor ---")
med_se = np.median([r["_se"] for r in ok])
z_bon = stats.norm.isf(0.5 * 0.05 / n)
z_80 = stats.norm.isf(0.20)
print("median SE(beta)      : %.4f  (per SD methylation, log-OR scale)" % med_se)
print("80%% power at Bonferroni alpha: beta >= %.3f  (OR %.2f)"
      % ((z_bon + z_80) * med_se, math.exp((z_bon + z_80) * med_se)))
print("observed |beta| range: %.3f - %.3f (median %.3f)"
      % (min(abs(r["_b"]) for r in ok), max(abs(r["_b"]) for r in ok),
         np.median([abs(r["_b"]) for r in ok])))
