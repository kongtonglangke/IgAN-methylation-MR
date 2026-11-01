"""Stage 1 gate test: cis-mQTL of candidate CpGs -> outcome (B-cell trait or IgAN).

INSTRUMENT SELECTION (the single most important design choice here)
------------------------------------------------------------------
GoDMC returns every significant cis-SNP per CpG, but those SNPs sit in strong
mutual LD (median 113 per CpG, up to 2322). Two failure modes follow:

  1. Feeding all of them to IVW treats correlated SNPs as independent
     observations. That shrinks the SE by ~sqrt(n_effective) and yields
     spurious p values around 1e-30. Observed live, see logs.
  2. GoDMC's own `clumped == True` flag cannot be used as the instrument filter
     either: it marks region-level LD index SNPs, which are frequently *not*
     the significant ones. For cg01140143 the strongest cis-SNP (p = 3.3e-18)
     is flagged clumped == False, while the 3 flagged SNPs are all non-
     significant.

We therefore use exactly ONE instrument per CpG: the strongest cis-SNP
(lowest p), analysed by Wald ratio. This is the standard drug-target MR
convention when cis instruments are few, it is immune to LD inflation, and it
keeps F-statistics high (median lead F ~ 195 in this dataset). Loss of
efficiency versus a multi-SNP score is small because the top SNP carries most
of the mQTL R2.

MATCHING NOTE
-------------
Orru 2020 GWAS-SSF files have an empty variant_id column (all NA), so SNP
identity must be matched on (chromosome, base_pair_location) in GRCh37, not
on rsID. GoDMC reports positions as 'chr7:7095466:SNP'.

Usage: python mr_cpg_bcell.py <trait_label> [max_cpg] [max_dist_kb]
"""
import csv, gzip, json, math, os, statistics, sys
from collections import defaultdict

PVAL = 5e-8
FTHR = 10.0


def norm_chrom(s):
    s = str(s).strip()
    if s.lower().startswith("chr"):
        s = s[3:]
    return s.upper()


def parse_snp(snp):
    """'chr7:7095466:SNP' -> ('7', 7095466)"""
    p = str(snp).split(":")
    if len(p) < 2:
        return None
    try:
        return (norm_chrom(p[0]), int(p[1]))
    except ValueError:
        return None


def load_candidates(max_cpg=None, max_dist=None):
    rows = []
    with open("results/cpg_instrument_strength.tsv") as fh:
        for r in csv.DictReader(fh, delimiter=chr(9)):
            if not r["status"].startswith("ok"):
                continue
            if float(r["lead_F"]) < FTHR:
                continue
            if max_dist is not None and int(r["dist_to_gene"]) > max_dist:
                continue
            rows.append(r)
    rows.sort(key=lambda x: -float(x["lead_F"]))
    if max_cpg:
        rows = rows[:max_cpg]
    return rows


def gwas_header(path):
    with gzip.open(path, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            return line.rstrip().split(chr(9))
    return []


def idx(header, *names):
    low = [h.lower() for h in header]
    for n in names:
        if n.lower() in low:
            return low.index(n.lower())
    return None


def stream_extract(path, wanted):
    """Single pass over the outcome GWAS, keeping only `wanted` (chr,pos) keys."""
    hdr = gwas_header(path)
    ci = {
        "chr": idx(hdr, "chromosome", "chr"),
        "pos": idx(hdr, "base_pair_location", "pos", "position", "bp"),
        "a1": idx(hdr, "effect_allele", "a1", "alt"),
        "a2": idx(hdr, "other_allele", "a2", "ref"),
        "beta": idx(hdr, "beta", "b", "effect"),
        "se": idx(hdr, "standard_error", "se", "stderr"),
        "pval": idx(hdr, "p_value", "pval", "p"),
        "eaf": idx(hdr, "effect_allele_frequency", "eaf", "freq"),
    }
    miss = [k for k, v in ci.items() if v is None and k != "eaf"]
    if miss:
        print("  WARN missing cols:", miss, "header:", hdr[:14], flush=True)
        return {}
    out = {}
    maxi = max(v for v in ci.values() if v is not None)
    n = 0
    with gzip.open(path, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip().split(chr(9))
            if len(p) <= maxi:
                continue
            n += 1
            try:
                key = (norm_chrom(p[ci["chr"]]), int(p[ci["pos"]]))
            except ValueError:
                continue
            if key not in wanted:
                continue
            try:
                eaf = float(p[ci["eaf"]]) if ci["eaf"] is not None else None
            except ValueError:
                eaf = None
            try:
                out[key] = (p[ci["a1"]].upper(), p[ci["a2"]].upper(),
                            float(p[ci["beta"]]), float(p[ci["se"]]),
                            float(p[ci["pval"]]), eaf)
            except ValueError:
                continue
    print("  scanned=%d matched=%d/%d" % (n, len(out), len(wanted)), flush=True)
    return out


COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def harmonise(bx, ax1, ax2, by, ay1, ay2):
    """Return `by` re-signed so that it refers to the exposure's effect allele ax1."""
    if ax1 == ay1 and ax2 == ay2:
        return by
    if ax1 == ay2 and ax2 == ay1:
        return -by
    if COMP.get(ax1) == ay1 and COMP.get(ax2) == ay2:
        return by
    if COMP.get(ax1) == ay2 and COMP.get(ax2) == ay1:
        return -by
    return None


def is_palindromic(a1, a2):
    return COMP.get(a1) == a2


def two_sided_p(z):
    return math.erfc(abs(z) / math.sqrt(2))


def main():
    trait = sys.argv[1]
    max_cpg = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    max_dist_kb = int(sys.argv[3]) if len(sys.argv) > 3 else None
    gwas = "data/gwas/%s.tsv.gz" % trait
    cands = load_candidates(max_cpg=max_cpg,
                            max_dist=(max_dist_kb * 1000 if max_dist_kb else None))
    print("trait=%s candidates=%d" % (trait, len(cands)), flush=True)

    gene_cache = {}
    cpg_inst = {}
    wanted = set()
    for c in cands:
        gene = c["gene"]
        if gene not in gene_cache:
            fp = os.path.join("data/mqtl", gene + ".json")
            if not os.path.exists(fp):
                gene_cache[gene] = None
            else:
                with open(fp) as fh:
                    data = json.load(fh)
                by_cpg = defaultdict(list)
                for r in data:
                    if r.get("cistrans") is not True:
                        continue
                    by_cpg[r.get("cpg")].append(r)
                gene_cache[gene] = by_cpg
        by_cpg = gene_cache[gene]
        if by_cpg is None:
            continue
        lst = by_cpg.get(c["cpg"], [])
        best = None
        for r in lst:
            try:
                p = float(r["pval_are"]); b = float(r["beta_are_a1"]); s = float(r["se_are"])
            except (TypeError, ValueError, KeyError):
                continue
            if s <= 0 or p >= PVAL or b == 0:
                continue
            k = parse_snp(r.get("snp"))
            if k is None:
                continue
            if best is None or p < best[5]:
                try:
                    fq = float(r.get("freq_a1"))
                except (TypeError, ValueError):
                    fq = 0.5
                best = (k, str(r["a1"]).upper(), str(r["a2"]).upper(), b, s, p, fq)
        if best is None:
            continue
        cpg_inst[c["cpg"]] = (c["gene"], float(c["lead_F"]), best)
        wanted.add(best[0])

    print("  cpgs_with_instruments=%d unique_positions=%d" % (len(cpg_inst), len(wanted)), flush=True)
    if not wanted:
        print("  no instruments -- abort", flush=True)
        return

    outd = stream_extract(gwas, wanted)

    rows = []
    for cpg, (gene, leadF, ins) in cpg_inst.items():
        k, a1, a2, b, s, p, fq = ins
        tag = "%s:%d" % k
        if k not in outd:
            rows.append((cpg, gene, trait, tag, "", "", "", "", "%.1f" % leadF, "no_snp_in_outcome"))
            continue
        oa1, oa2, ob, os_, op, oeaf = outd[k]
        if is_palindromic(a1, a2) and oeaf is not None and 0.35 < oeaf < 0.65:
            rows.append((cpg, gene, trait, tag, "", "", "", "", "%.1f" % leadF, "palindromic_ambig"))
            continue
        hb = harmonise(b, a1, a2, ob, oa1, oa2)
        if hb is None:
            rows.append((cpg, gene, trait, tag, "", "", "", "", "%.1f" % leadF, "allele_mismatch"))
            continue
        if os_ <= 0:
            rows.append((cpg, gene, trait, tag, "", "", "", "", "%.1f" % leadF, "bad_se"))
            continue
        bhat = hb / b
        se = abs(os_ / b)
        z = bhat / se if se > 0 else 0.0
        rows.append((cpg, gene, trait, tag, "%.4f" % b, "%.4f" % bhat, "%.4f" % se,
                     "%.3e" % two_sided_p(z), "%.1f" % leadF, "ok"))

    os.makedirs("results", exist_ok=True)
    out = "results/mr_stage1_%s.tsv" % trait
    with open(out, "w") as f:
        f.write(chr(9).join(["cpg", "gene", "trait", "snp", "bx", "b", "se", "pval",
                             "lead_F", "status"]) + chr(10))
        for r in rows:
            f.write(chr(9).join(str(x) for x in r) + chr(10))
    ok = [r for r in rows if r[-1] == "ok"]
    sig = [r for r in ok if float(r[7]) < 0.05]
    print("  RESULT total=%d ok=%d p<0.05=%d (%.1f%%)" %
          (len(rows), len(ok), len(sig), 100.0 * len(sig) / max(1, len(ok))), flush=True)
    if ok:
        ps = sorted(float(r[7]) for r in ok)
        n = len(ok)
        D = max(max(abs((i + 1) / n - x), abs(x - i / n)) for i, x in enumerate(ps))
        lam = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * D
        ks = 2 * sum((-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam) for k in range(1, 60))
        print("  min_p=%.3e  median_p=%.3f" % (ps[0], statistics.median(ps)), flush=True)
        print("  KS vs U(0,1): D=%.4f p=%.3f  (p>0.05 => indistinguishable from null)" %
              (D, min(1, max(0, ks))), flush=True)
        print("  expected p<0.05 under null = %.1f, observed = %d" %
              (0.05 * n, sum(1 for x in ps if x < 0.05)), flush=True)


if __name__ == "__main__":
    main()
