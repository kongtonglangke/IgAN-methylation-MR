"""Fetch IgAN (Sakaue 2021, GCST90018866) summary statistics for the cis-mQTL
instrument SNPs, using tabix random access instead of a 540 MB bulk download.

Writes data/gwas/igan_lookup.json : {"chr:pos": [field_list, ...]}
"""
import sys, csv, json, time

sys.path.insert(0, "scripts")
from tabix_http import TabixHTTP

BASE = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90018001-GCST90019000/GCST90018866/harmonised/")

rows = [r for r in csv.DictReader(open("results/mr_stage1_bcell_igdp.tsv"), delimiter=chr(9))]
pos = set()
for r in rows:
    if ":" in r["snp"]:
        c, p = r["snp"].split(":")[:2]
        pos.add((c, p))
print("positions to fetch: %d" % len(pos), flush=True)

tb = TabixHTTP(BASE + "GCST90018866.h.tsv.gz", BASE + "GCST90018866.h.tsv.gz.tbi")
t0 = time.time()
res = tb.fetch(pos, nthreads=8)
tot = len(res)
hit = sum(1 for k, v in res.items() if isinstance(v, list) and v)
err = [k for k, v in res.items() if not isinstance(v, list)]
print("fetched %.1fs  total=%d hit=%d (%.1f%%) err=%d"
      % (time.time() - t0, tot, hit, 100.0 * hit / max(1, tot), len(err)), flush=True)

out = {("%s:%s" % k): v for k, v in res.items() if isinstance(v, list)}
json.dump(out, open("data/gwas/igan_lookup.json", "w"))
print("saved data/gwas/igan_lookup.json", flush=True)
