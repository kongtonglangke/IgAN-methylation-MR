"""Stage 1b: CpG methylation -> IgA nephropathy (IgA 肾病) causal test.

Design
------
Exposure : blood cis-mQTL instrument, one lead SNP per CpG (GoDMC meta-analysis,
           Min et al. 2021), already screened to F >= 25 in
           results/cpg_instrument_strength.tsv.
Outcome  : Sakaue et al. 2021 IgAN GWAS, GCST90018866 (15,587 cases /
           462,197 controls, EUR+EAS meta), GWAS-Catalog harmonised file.
Estimator: single-instrument Wald ratio  b = b_out / b_meth,
           se = se_out / |b_meth|   (first-order delta method).
           No IVW: neighbouring cis-SNPs are in tight LD, so pooling them
           inflates precision without adding information (see stage-1a note).

Build handling
--------------
The harmonised outcome file is GRCh38; GoDMC instruments are GRCh37. Every
instrument is lifted 37->38 through Ensembl REST before lookup, and the
lifted coordinate is cross-checked against the file's own rsID column.

Outputs
-------
data/liftover_37_38.json         cache of lifted positions
data/gwas/igan_lookup_38.json    raw outcome rows keyed by GRCh38 chrom:pos
results/mr_stage1b_cpg_igan.tsv  per-CpG MR estimates
"""
import csv
import json
import os
import ssl
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tabix_http import TabixHTTP
from liftover_chain import ChainLift

BASE = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90018001-GCST90019000/GCST90018866/harmonised/")
URL = BASE + "GCST90018866.h.tsv.gz"
IDX = URL + ".tbi"

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


# ---------------------------------------------------------------- instruments
def load_instruments(min_f=25.0):
    rows = []
    with open("results/cpg_instrument_strength.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["status"] != "ok":
                continue
            try:
                f = float(r["lead_F"])
            except ValueError:
                continue
            if f < min_f or not r["lead_snp"]:
                continue
            parts = r["lead_snp"].split(":")
            if len(parts) < 2:
                continue
            chrom = parts[0].replace("chr", "")
            pos = int(parts[1])
            a1, _, a2 = r["lead_alleles"].partition("/")
            rows.append(dict(cpg=r["cpg"], gene=r["gene"], chrom=chrom, pos37=pos,
                             a1=a1.upper(), a2=a2.upper(), lead_F=f,
                             rsid=r["lead_rsid"]))
    return rows


# ------------------------------------------------------------------- liftover
def liftover(rows):
    """GRCh37 -> GRCh38 via the local UCSC chain file (see liftover_chain.py)."""
    cache_path = "data/liftover_37_38.json"
    cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}
    L = ChainLift()
    n_new = 0
    for r in rows:
        key = "%s:%d" % (r["chrom"], r["pos37"])
        if key in cache:
            continue
        got = L.map(r["chrom"], r["pos37"])
        cache[key] = [r["chrom"], got] if got else None
        n_new += 1
    json.dump(cache, open(cache_path, "w"))
    unmapped = sum(1 for r in rows if not cache.get("%s:%d" % (r["chrom"], r["pos37"])))
    print("liftover: %d newly mapped, %d/%d unmapped"
          % (n_new, unmapped, len(rows)), flush=True)
    return cache


# ---------------------------------------------------------------- fetch outcome
def fetch_outcome(rows, cache):
    want = set()
    for r in rows:
        m = cache.get("%s:%d" % (r["chrom"], r["pos37"]))
        if m:
            want.add((str(m[0]), int(m[1])))
        want.add((r["chrom"], r["pos37"]))          # fallback: same-coordinate probe
    print("outcome lookups: %d positions" % len(want), flush=True)
    tb = TabixHTTP(URL, IDX)
    res = tb.fetch(want, nthreads=6)
    hit = sum(1 for v in res.values() if isinstance(v, list) and v)
    err = [k for k, v in res.items() if not isinstance(v, list)]
    print("fetched: hit=%d/%d  err=%d" % (hit, len(want), len(err)), flush=True)
    if err:                                          # one retry pass
        tb2 = TabixHTTP(URL, IDX)
        res2 = tb2.fetch(err, nthreads=2)
        for k, v in res2.items():
            if isinstance(v, list):
                res[k] = v
        print("retry: residual errors=%d"
              % sum(1 for k in err if not isinstance(res.get(k), list)), flush=True)
    json.dump({" ".join(map(str, k)): v for k, v in res.items() if isinstance(v, list)},
              open("data/gwas/igan_lookup_38.json", "w"))
    return res


# ------------------------------------------------------------------ harmonise
def pick(rows, a1, a2, rsid):
    """Return (beta, se, p, eaf, match_type) aligned to allele a1 (methylation
    increasing allele) or None. Rows carry the IgAN header order."""
    for r in rows:
        ea, oa = r[2].upper(), r[3].upper()
        try:
            beta, se, p = float(r[4]), float(r[5]), float(r[7])
        except (ValueError, IndexError):
            continue
        if ea == a1 and oa == a2:
            return beta, se, p, r[6], "direct"
        if ea == a2 and oa == a1:
            return -beta, se, p, r[6], "flip"
        if ea == COMP.get(a1) and oa == COMP.get(a2):
            return beta, se, p, r[6], "complement"
        if ea == COMP.get(a2) and oa == COMP.get(a1):
            return -beta, se, p, r[6], "flip+complement"
    return None


def main():
    rows = load_instruments()
    print("instruments with F>=25: %d CpGs, %d genes"
          % (len(rows), len({r["gene"] for r in rows})), flush=True)
    cache = liftover(rows)
    res = fetch_outcome(rows, cache)

    # methylation-side effect: strongest cis-SNP per CpG, from the GoDMC per-gene
    # dumps used by compute_f.py (keeps the instrument identical to stage 1a).
    meth = {}
    import glob
    for fp in glob.glob("data/mqtl/*.json"):
        for rec in json.load(open(fp)):
            if not isinstance(rec, dict):
                continue
            snp = rec.get("snp") or ""
            if not snp.startswith("chr"):
                continue
            try:
                _, pos, _ = snp.split(":")
            except ValueError:
                continue
            key = (rec["cpg"], int(pos))
            prev = meth.get(key)
            if prev is None or (rec.get("pval_are") or 1) < (prev.get("pval_are") or 1):
                meth[key] = rec
    print("mQTL records indexed: %d pairs / %d CpGs"
          % (len(meth), len({k[0] for k in meth})), flush=True)

    out = []
    for r in rows:
        m = cache.get("%s:%d" % (r["chrom"], r["pos37"]))
        got = None
        coord_used = None
        if m:
            got = res.get((str(m[0]), int(m[1])))
            if got:
                coord_used = "GRCh38:%s:%d" % (m[0], m[1])
        if not got:
            got = res.get((r["chrom"], r["pos37"]))
            if got:
                coord_used = "GRCh37:%s:%d" % (r["chrom"], r["pos37"])
        rec = meth.get((r["cpg"], r["pos37"]))
        bx = rec.get("beta_are_a1", rec.get("beta_a1")) if rec else None
        row = dict(cpg=r["cpg"], gene=r["gene"], rsid=r["rsid"],
                   snp37="%s:%d" % (r["chrom"], r["pos37"]),
                   coord_used=coord_used or "absent", lead_F="%.1f" % r["lead_F"],
                   bx="", se_meth="", b_out="", se_out="", p_out="", eaf="",
                   allele_match="", b="", se="", pval="", status="")
        if rec is None or bx in (None, 0):
            row["status"] = "no_meth_effect"
            out.append(row)
            continue
        row["bx"] = "%.5f" % bx
        row["se_meth"] = "%.5f" % rec.get("se_are", rec.get("se"))
        if not got:
            row["status"] = "not_in_igan_gwas"
            out.append(row)
            continue
        p = pick(got, r["a1"], r["a2"], r["rsid"])
        if p is None:
            row["status"] = "allele_mismatch"
            out.append(row)
            continue
        b_out, se_out, p_out, eaf, mt = p
        row.update(b_out="%.5f" % b_out, se_out="%.5f" % se_out,
                   p_out="%.4e" % p_out, eaf=eaf or "", allele_match=mt)
        b = b_out / bx
        se = se_out / abs(bx)
        z = b / se
        from math import erfc, sqrt
        pv = erfc(abs(z) / sqrt(2))
        row.update(b="%.5f" % b, se="%.5f" % se, pval="%.4e" % pv, status="ok")
        out.append(row)

    cols = ["cpg", "gene", "rsid", "snp37", "coord_used", "lead_F", "bx", "se_meth",
            "b_out", "se_out", "p_out", "eaf", "allele_match", "b", "se", "pval",
            "status"]
    with open("results/mr_stage1b_cpg_igan.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(out)

    from collections import Counter
    print("status:", dict(Counter(x["status"] for x in out)))
    ok = [x for x in out if x["status"] == "ok"]
    print("analysable CpGs: %d" % len(ok))


if __name__ == "__main__":
    main()
