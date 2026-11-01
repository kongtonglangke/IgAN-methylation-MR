"""Is the DEFA1 mQTL a copy-number signature or a local regulatory effect?

A copy-number variant that the SNP tags would shift methylation of essentially
every CpG across the duplicated block, in the same direction, with the same
lead variant.  A local regulatory effect hits a subset of CpGs and different
CpGs usually nominate different lead SNPs.  Both signatures are measurable
from the GoDMC per-gene dumps we already hold.

Controls: the HORMAD2 and CFH loci, which are not segmental duplications.
"""
import csv
import json
import os
import statistics
from collections import defaultdict

GENES = {
    "DEFA1":    ("rs2738048", "chr8"),
    "HORMAD2":  ("rs12537",   "chr22"),
    "CFH":      ("rs6677604", "chr1"),
}

SIG = 5e-8
SUB = 1e-3


def load(gene):
    fp = "data/mqtl/%s.json" % gene
    if not os.path.exists(fp):
        return None
    return json.load(open(fp))


def snp_of(rec):
    s = rec.get("snp") or ""
    parts = s.split(":")
    return parts[1] if len(parts) == 3 else None


def summarise(gene, lead_rsid):
    recs = load(gene)
    if recs is None:
        print("%-9s  (no local dump)" % gene)
        return None
    by_cpg = defaultdict(list)
    for r in recs:
        by_cpg[r["cpg"]].append(r)
    n_cpg = len(by_cpg)

    # ---- per-CpG strongest instrument, ignoring the lead SNP
    own_lead = {}
    for cpg, rs in by_cpg.items():
        best = None
        for r in rs:
            p = r.get("pval_mre")
            if p is None:
                continue
            if best is None or p < best[0]:
                best = (p, r.get("rsid") or snp_of(r))
        if best:
            own_lead[cpg] = best

    # ---- how the lead IgAN SNP behaves across all CpGs
    hits = []
    for cpg, rs in by_cpg.items():
        for r in rs:
            if r.get("rsid") == lead_rsid:
                hits.append((cpg, r))
                break

    print("=" * 74)
    print("%s   lead SNP %s   CpGs in window: %d" % (gene, lead_rsid, n_cpg))
    print("=" * 74)
    if not hits:
        print("  lead SNP absent from this window's dump")
        print("  CpGs with a resolvable own-lead SNP: %d" % len(own_lead))
        return None

    sig = [(c, r) for c, r in hits if (r.get("pval_mre") or 1) < SIG]
    sub = [(c, r) for c, r in hits if (r.get("pval_mre") or 1) < SUB]
    betas = [r["beta_a1"] for c, r in sig]
    isq = [r.get("hetisq") or 0 for c, r in hits]
    pos = [(int(snp_of(r) or 0), r["beta_a1"]) for c, r in sig]

    print("  CpGs carrying a record for this SNP      : %d / %d" % (len(hits), n_cpg))
    print("  ... significant at p<5e-8                : %d (%.0f%%)"
          % (len(sig), 100.0 * len(sig) / max(len(hits), 1)))
    print("  ... nominal at p<1e-3                    : %d (%.0f%%)"
          % (len(sub), 100.0 * len(sub) / max(len(hits), 1)))
    if betas:
        pos_b = sum(1 for b in betas if b > 0)
        print("  effect direction among significant       : %d up / %d down"
              % (pos_b, len(betas) - pos_b))
        print("  |beta| median %.3f  range %.3f - %.3f"
              % (statistics.median(abs(b) for b in betas),
                 min(abs(b) for b in betas), max(abs(b) for b in betas)))
    print("  heterogeneity I2 median %.1f  max %.1f  (n with I2>75: %d)"
          % (statistics.median(isq), max(isq), sum(1 for x in isq if x > 75)))

    # ---- do the CpGs share ONE lead SNP, or many?
    lead_ids = [v[1] for v in own_lead.values() if v[1]]
    cnt = defaultdict(int)
    for l in lead_ids:
        cnt[l] += 1
    top = sorted(cnt.items(), key=lambda kv: -kv[1])[:5]
    print("\n  distinct own-lead SNPs across CpGs       : %d for %d CpGs"
          % (len(cnt), len(lead_ids)))
    print("  most frequent own-lead SNPs:")
    for k, v in top:
        tag = "  <-- IgAN lead" if k == lead_rsid else ""
        print("     %-14s %3d CpGs%s" % (k, v, tag))

    # ---- position span of the affected CpGs
    if pos:
        ps = sorted(pos)
        print("\n  affected CpGs span %s - %s (by lead-SNP position)"
              % (ps[0][0], ps[-1][0]))
        # contiguity: are the affected CpGs clustered or scattered?
        print("  (a CNV would shift every CpG across the duplicated block;")
        print("   a local effect concentrates on a subset)")
    return dict(gene=gene, lead=lead_rsid, n_cpg=n_cpg, n_hits=len(hits),
                n_sig=len(sig), n_sub=len(sub),
                med_abs_beta=statistics.median(abs(b) for b in betas) if betas else None,
                n_distinct_lead=len(cnt))


def main():
    out = []
    for gene, (lead, _) in GENES.items():
        r = summarise(gene, lead)
        if r:
            out.append(r)
        print()
    if out:
        json.dump(out, open("results/defa1_cnv_signature.json", "w"), indent=1)
        print("saved results/defa1_cnv_signature.json")


if __name__ == "__main__":
    main()
