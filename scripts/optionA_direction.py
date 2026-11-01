#!/usr/bin/env python
"""Option A, step 3: align the IgAN risk allele with the GoDMC effect allele and
report the direction in which each risk allele shifts CpG methylation.

The GoDMC beta (beta_a1) is per copy of allele a1; the GWAS Catalog OR is per
copy of the reported effect allele.  If the two alleles differ the sign of the
mQTL effect must be flipped before it can be interpreted as "the risk allele
raises/lowers methylation".
"""
import csv
import json
import os
from collections import defaultdict

CANON = {"A", "C", "G", "T"}


def load_mqtl(gene, cache={}):
    if gene not in cache:
        path = "data/mqtl/%s.json" % gene
        cache[gene] = json.load(open(path)) if os.path.exists(path) else None
    return cache[gene]


def norm(a):
    return (a or "").strip().upper()


def main():
    # effect allele reported for each IgAN variant; a rsID can appear in several
    # association records (different studies), so a missing allele must never
    # shadow a populated one
    eff = {}
    for r in csv.DictReader(open("data/gwas/igan_catalog_assoc.tsv"), delimiter="\t"):
        rs, ea = r["rsid"], norm(r["effect_allele"])
        if rs and ea in CANON and not eff.get(rs):
            eff[rs] = ea

    pairs = list(csv.DictReader(open("results/optionA_igan_variant_mqtl.tsv"),
                                delimiter="\t"))
    out = []
    for r in pairs:
        recs = load_mqtl(r["gene"])
        if not recs:
            continue
        hits = [x for x in recs if x.get("rsid") == r["igan_rsid"] and x["cpg"] == r["cpg"]]
        if not hits:
            continue
        x = hits[0]
        p_re = x.get("pval_mre")
        if p_re is None or p_re >= 5e-8:
            continue
        a1, a2 = norm(x.get("a1")), norm(x.get("a2"))
        ea = eff.get(r["igan_rsid"], "")
        if ea and ea in CANON and a1 in CANON and a2 in CANON:
            if ea == a1:
                align, sign = "effect_allele_is_a1", 1
            elif ea == a2:
                align, sign = "effect_allele_is_a2", -1
            elif a1 == a2:
                align, sign = "ambiguous_palindrome", 0
            else:
                align, sign = "allele_mismatch", 0
        else:
            align, sign = "unresolved", 0
        beta = x["beta_a1"] * sign if sign else None
        out.append({
            "gene": r["gene"], "rsid": r["igan_rsid"], "cpg": r["cpg"],
            "igan_p": r["igan_p"], "igan_or": r["igan_or"],
            "effect_allele": ea, "mqtl_a1": a1, "mqtl_a2": a2,
            "alignment": align,
            "beta_per_risk_allele": "" if beta is None else round(beta, 4),
            "direction": ("" if beta is None else
                          ("risk allele METHYLATES" if beta > 0 else "risk allele DEMETHYLATES")),
            "p_random_effects": p_re, "het_isq": x.get("hetisq"),
            "n_studies": x.get("num_studies"),
        })

    out.sort(key=lambda r: (r["gene"], r["rsid"], float(r["p_random_effects"])))
    os.makedirs("results", exist_ok=True)
    with open("results/optionA_direction.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    print("significant (random-effects) variant-CpG pairs with resolvable alleles: %d"
          % len(out))
    from collections import Counter
    print("alignment:", dict(Counter(r["alignment"] for r in out)))
    print()

    # one representative CpG per (variant, gene): strongest random-effects p
    best = {}
    for r in out:
        k = (r["rsid"], r["gene"])
        if k not in best or float(r["p_random_effects"]) < float(best[k]["p_random_effects"]):
            best[k] = r
    print("%-14s %-10s %-11s %-6s %-6s %9s %-9s %s" % (
        "rsid", "gene", "top_cpg", "EA", "beta", "p_RE", "I2", "direction"))
    for k, r in sorted(best.items(), key=lambda kv: kv[0][1]):
        print("%-14s %-10s %-11s %-6s %-6s %9.1e %-9s %s" % (
            r["rsid"], r["gene"], r["cpg"], r["effect_allele"] or "?",
            r["beta_per_risk_allele"], float(r["p_random_effects"]),
            r["het_isq"], r["direction"] or r["alignment"]))

    up = sum(1 for r in best.values() if r["direction"].endswith("METHYLATES"))
    dn = sum(1 for r in best.values() if r["direction"].endswith("DEMETHYLATES"))
    print("\ndirection tally (representative): %d methylating, %d demethylating" % (up, dn))


if __name__ == "__main__":
    main()
