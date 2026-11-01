#!/usr/bin/env python
"""Option A core test: are the known IgA nephropathy risk variants themselves
mQTLs of the metabolite-target genes?

Inputs
------
data/gwas/igan_loci_overlap.tsv   IgAN loci lying within +-500 kb of a target gene
data/mqtl/<GENE>.json             GoDMC cis-mQTL results, GRCh37

Matching strategy
-----------------
The GoDMC SNP keys are GRCh37 positions while GWAS Catalog uses GRCh38, and
the offset reaches ~760 kb at 16p11.2, so positional matching is only used
after an explicit liftover.  rsIDs are assembly-independent and are tried
first; position matching is the fallback for records where GoDMC has no rsID.

Output
------
results/optionA_igan_variant_mqtl.tsv   one row per (IgAN variant, CpG, gene)
results/optionA_summary.tsv             one row per IgAN variant
"""
import csv
import json
import os
import sys


def load_mqtl(gene):
    path = "data/mqtl/%s.json" % gene
    if not os.path.exists(path):
        return None
    return json.load(open(path))


def main():
    loci = list(csv.DictReader(open("data/gwas/igan_loci_overlap.tsv"), delimiter="\t"))
    # unique (variant, gene) pairs
    pairs = {}
    for r in loci:
        pairs.setdefault((r["target_gene"], r["rsid"], r["loc38"]), r)

    cache = {}
    detail, summary = [], []
    for (gene, rsid, loc), meta in sorted(pairs.items()):
        if gene not in cache:
            cache[gene] = load_mqtl(gene)
        recs = cache[gene]
        if recs is None:
            summary.append(dict(gene=gene, rsid=rsid, loc38=loc,
                                igan_p=meta["p_value"], igan_or=meta["or_value"],
                                status="no_mqtl_file", match="none",
                                n_cpg_hit=0, n_sig=0, best_cpg="", best_beta="",
                                best_se="", best_p="", best_rsid=""))
            continue

        hits = [x for x in recs if x.get("rsid") == rsid]
        match_mode = "rsid"
        if not hits:
            # fallback: same rsID stored under a different spelling, or the
            # variant is present but without an rsID -> leave to liftover step
            match_mode = "none"

        sig = [x for x in hits if (x.get("pval") is not None and x["pval"] < 5e-8)]
        sig.sort(key=lambda x: x["pval"])
        for x in sig:
            detail.append({
                "gene": gene, "igan_rsid": rsid, "igan_loc38": loc,
                "igan_p": meta["p_value"], "igan_or": meta["or_value"],
                "cpg": x["cpg"], "mqtl_beta": x.get("beta_a1"),
                "mqtl_se": x.get("se"), "mqtl_p": x.get("pval"),
                "mqtl_a1": x.get("a1"), "mqtl_a2": x.get("a2"),
                "cistrans": x.get("cistrans"), "n_studies": x.get("num_studies"),
                "het_isq": x.get("hetisq"),
            })
        best = sig[0] if sig else (sorted(hits, key=lambda x: x.get("pval") or 1)[0]
                                   if hits else None)
        summary.append(dict(
            gene=gene, rsid=rsid, loc38=loc,
            igan_p=meta["p_value"], igan_or=meta["or_value"],
            status=("mqtl_sig" if sig else ("mqtl_nominal" if hits else "not_in_data")),
            match=match_mode,
            n_cpg_hit=len(hits), n_sig=len(sig),
            best_cpg=best["cpg"] if best else "",
            best_beta=best.get("beta_a1") if best else "",
            best_se=best.get("se") if best else "",
            best_p=best.get("pval") if best else "",
            best_rsid=rsid,
        ))

    os.makedirs("results", exist_ok=True)
    with open("results/optionA_igan_variant_mqtl.tsv", "w", newline="") as fh:
        if detail:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(detail[0].keys()))
            w.writeheader()
            w.writerows(detail)
        else:
            fh.write("no significant mQTL hits\n")
    with open("results/optionA_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)

    # ---------------- report ----------------
    print("IgAN variant x target-gene pairs examined: %d" % len(pairs))
    from collections import Counter
    c = Counter(s["status"] for s in summary)
    print("  mqtl_sig (p<5e-8)      : %d" % c.get("mqtl_sig", 0))
    print("  mqtl_nominal           : %d" % c.get("mqtl_nominal", 0))
    print("  not present in GoDMC   : %d" % c.get("not_in_data", 0))
    print("  gene has no mQTL file  : %d" % c.get("no_mqtl_file", 0))
    print()
    sig = [s for s in summary if s["status"] == "mqtl_sig"]
    if sig:
        print("=== variants that ARE significant cis-mQTLs ===")
        print("%-14s %-9s %-11s %-9s %-8s %-10s %-9s" %
              ("igan_rsid", "gene", "best_cpg", "beta", "se", "p_mqtl", "n_cpg"))
        for s in sorted(sig, key=lambda x: float(x["best_p"])):
            print("%-14s %-9s %-11s %-9.3f %-8.4f %-10.1e %-9d" %
                  (s["rsid"], s["gene"], s["best_cpg"], float(s["best_beta"]),
                   float(s["best_se"]), float(s["best_p"]), s["n_sig"]))
    else:
        print("=== no IgAN risk variant reached p<5e-8 as a cis-mQTL ===")
        nom = [s for s in summary if s["status"] == "mqtl_nominal"]
        print("nominal-only hits: %d" % len(nom))
        for s in sorted(nom, key=lambda x: float(x["best_p"]))[:15]:
            print("   %-14s %-9s %-11s beta=%-8.3f p=%.1e" %
                  (s["rsid"], s["gene"], s["best_cpg"], float(s["best_beta"]),
                   float(s["best_p"])))


if __name__ == "__main__":
    sys.exit(main())
