"""Diagnose whether GCST90018866 (Sakaue 2021 IgAN) harmonised file is
genome-wide or filtered to significant variants.

Counts rows retrieved from a few 1-Mb windows and prints the p-value
distribution; a filtered file shows only a handful of near-significant rows
per window plus a huge void elsewhere.
"""
import sys, time

sys.path.insert(0, "scripts")
from tabix_http import TabixHTTP

BASE = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90018001-GCST90019000/GCST90018866/harmonised/")
URL = BASE + "GCST90018866.h.tsv.gz"
IDX = URL + ".tbi"

tb = TabixHTTP(URL, IDX)
print("index loaded: %d refs" % len(tb.refs), flush=True)

wins = [("1", 196600000, 197700000),   # 1q32 CFH locus, should be dense
        ("22", 21800000, 22900000),    # 22q12 HORMAD2 locus
        ("6", 25000000, 35000000),     # MHC
        ("9", 5000000, 6000000),       # arbitrary background window
        ("4", 100000000, 101000000)]   # arbitrary background window
for c, s, e in wins:
    t0 = time.time()
    rows = tb.query(c, s, e, max_blocks=4000)
    pv = []
    for f in rows:
        try:
            pv.append(float(f[7]))
        except Exception:
            pass
    n_sig = sum(1 for p in pv if p < 5e-8)
    n_sug = sum(1 for p in pv if p < 1e-5)
    print("%s:%d-%d  rows=%-6d  p<5e-8=%-4d p<1e-5=%-4d  min_p=%s  (%.1fs)"
          % (c, s, e, len(rows), n_sig, n_sug,
             ("%.2e" % min(pv)) if pv else "NA", time.time() - t0), flush=True)
