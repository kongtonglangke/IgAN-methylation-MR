"""Compute per-CpG instrument strength from GoDMC cis-mQTL.

Key corrections vs naive version:
 1. Only SNPs with cistrans == True are used (GoDMC cis definition = within +/-1Mb).
 2. CpGs whose SNPs are ALL trans are dropped -- they cannot support cis-MR.
 3. The GoDMC range endpoint can return CpGs outside the queried window, so the
    CpG position is estimated as the median position of its cis SNPs and the CpG
    is discarded if that estimate falls far outside the gene window.
 4. Lead SNP chosen among cis SNPs only; F = (beta_are_a1 / se_are)^2,
    ARE model per Min et al. 2021 Nat Genet.

Output: results/cpg_instrument_strength.tsv
"""
import glob, json, os, statistics
from collections import defaultdict

PVAL = 5e-8
FTHR = 10.0
CIS_TOL = 1500000   # allow 1.5 Mb slack outside gene window when judging membership
MHC_LO, MHC_HI = 25000000, 35000000

genes = {}
with open("data/gene_coords_grch37.tsv") as fh:
    next(fh)
    for line in fh:
        p = line.rstrip().split(chr(9))
        if len(p) >= 4 and p[1] != "NA":
            genes[p[0]] = (p[1], int(p[2]), int(p[3]))

WIN = 500000
rows_out = []
files = sorted(glob.glob("data/mqtl/*.json"))
print("files=%d" % len(files), flush=True)

for fp in files:
    gene = os.path.basename(fp).replace(".json", "")
    if gene not in genes:
        continue
    chrom, gstart, gend = genes[gene]
    lo = max(1, gstart - WIN - CIS_TOL)
    hi = gend + WIN + CIS_TOL
    with open(fp) as fh:
        data = json.load(fh)
    by_cpg = defaultdict(list)
    for r in data:
        c = r.get("cpg")
        if c:
            by_cpg[c].append(r)

    n_ok = n_transonly = n_out = 0
    for cpg, lst in by_cpg.items():
        tot = len(lst)
        cis = [r for r in lst if r.get("cistrans") is True]
        if not cis:
            rows_out.append((cpg, gene, tot, 0, 0, "", "", "", "", "", 0.0, 0, "trans_only"))
            n_transonly += 1
            continue
        pos = []
        for r in cis:
            try:
                pos.append(int(str(r.get("snp", "")).split(":")[1]))
            except (IndexError, ValueError):
                pass
        if not pos:
            rows_out.append((cpg, gene, tot, len(cis), 0, "", "", "", "", "", 0.0, 0, "no_pos"))
            continue
        cpg_est = int(statistics.median(pos))
        if cpg_est < lo or cpg_est > hi:
            rows_out.append((cpg, gene, tot, len(cis), 0, "", "", "", str(cpg_est), "", 0.0, 0, "out_of_window"))
            n_out += 1
            continue

        sig = []
        n_f10 = 0
        for r in cis:
            b = r.get("beta_are_a1")
            s = r.get("se_are")
            p = r.get("pval_are")
            try:
                b = float(b); s = float(s); p = float(p)
            except (TypeError, ValueError):
                continue
            if s <= 0:
                continue
            f = (b / s) ** 2
            if f > FTHR:
                n_f10 += 1
            if p < PVAL:
                sig.append((p, b, s, f, r))
        if not sig:
            rows_out.append((cpg, gene, tot, len(cis), 0, "", "", "", str(cpg_est), "", 0.0, n_f10, "no_sig_cis"))
            continue
        sig.sort(key=lambda x: x[0])
        p, b, s, f, best = sig[0]
        dist = 0
        if gstart <= cpg_est <= gend:
            dist = 0
        else:
            dist = min(abs(cpg_est - gstart), abs(cpg_est - gend))
        in_mhc = "MHC" if (chrom == "6" and MHC_LO <= cpg_est <= MHC_HI) else ""
        rows_out.append((
            cpg, gene, tot, len(cis), len(sig),
            best.get("rsid", ""), best.get("snp", ""),
            str(best.get("a1", "")) + "/" + str(best.get("a2", "")),
            str(cpg_est), str(dist),
            round(f, 2), n_f10, "ok" + ("|" + in_mhc if in_mhc else ""),
        ))
        n_ok += 1
    print("%-12s cpgs=%-4d ok=%-4d trans_only=%-3d out=%-3d" % (gene, len(by_cpg), n_ok, n_transonly, n_out), flush=True)

os.makedirs("results", exist_ok=True)
hdr = ["cpg", "gene", "n_snp_total", "n_cis", "n_cis_sig", "lead_rsid", "lead_snp",
       "lead_alleles", "cpg_pos_est", "dist_to_gene", "lead_F", "n_cis_F10", "status"]
with open("results/cpg_instrument_strength.tsv", "w") as f:
    f.write(chr(9).join(hdr) + chr(10))
    for r in rows_out:
        f.write(chr(9).join(str(x) for x in r) + chr(10))

n_tot = len(rows_out)
n_pass = sum(1 for r in rows_out if r[12].startswith("ok") and r[10] >= FTHR)
print("=" * 60)
print("cpg_total=%d  trans_only=%d  out_of_window=%d  usable_F>=10=%d"
      % (n_tot,
         sum(1 for r in rows_out if r[12] == "trans_only"),
         sum(1 for r in rows_out if r[12] == "out_of_window"),
         n_pass), flush=True)
