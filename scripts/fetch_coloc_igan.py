#!/usr/bin/env python
"""Fetch IgAN summary statistics around each colocalisation region.

CpG positions come from results/cpg_instrument_strength.tsv (cpg_pos_est,
GRCh37) with a fallback to the GoDMC /info endpoint cached in cpg_info.json.
Positions are lifted to GRCh38 with the UCSC chain so the tabix query lands on
the right coordinates in the Sakaue 2021 file (GCST90018866, GRCh38).

Because the CpGs of one region sit within a few hundred kb of each other, one
window per region is fetched instead of one per CpG — 7 queries instead of 19.

Output: data/coloc/igan_<region>.tsv
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from liftover_chain import ChainLift          # noqa: E402
from tabix_http import TabixHTTP             # noqa: E402

URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
       "GCST90018001-GCST90019000/GCST90018866/harmonised/GCST90018866.h.tsv.gz")
IDX = URL + ".tbi"
WINDOW = 750_000
OUT = "data/coloc"
FIELDS = ["rsid", "chr", "pos38", "ea", "oa", "beta", "se", "eaf", "p"]

REGIONS = {
    "CFH_CFHR_1q32":     ["cg13682187", "cg11348583", "cg13441112", "cg07209298"],
    "DEFA1_8p23":        ["cg00705004", "cg07297068", "cg19304352", "cg05969298"],
    "ITGAM_ITGAX_16p11": ["cg03100639", "cg10451425", "cg02846316", "cg04018474"],
    "HORMAD2_22q12":     ["cg27665648", "cg00539347"],
    "MADCAM1_19p13":     ["cg23386659", "cg27535410"],
    "ST6GAL1_3q27":      ["cg17448192", "cg25599673"],
    "TNFRSF13B_17p11":   ["cg02671941"],
}
CPGS = {c for v in REGIONS.values() for c in v}


def load_cpg_pos():
    pos = {}
    gene_chr = {}
    with open("data/gene_coords_grch37.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gene_chr[r["gene"]] = str(r["chr"]).replace("chr", "")
    for r in csv.DictReader(open("results/cpg_instrument_strength.tsv"),
                            delimiter="\t"):
        if r["cpg"] in CPGS and r.get("cpg_pos_est"):
            c = gene_chr.get(r["gene"])
            if c:
                pos[r["cpg"]] = (c, int(float(r["cpg_pos_est"])))

    info_path = os.path.join(OUT, "cpg_info.json")
    if os.path.exists(info_path):
        info = json.load(open(info_path))
        for cpg in CPGS:
            if cpg in pos or cpg not in info:
                continue
            d = info[cpg]
            d = d[0] if isinstance(d, list) else d
            if d.get("chr") and d.get("pos"):
                pos[cpg] = (str(d["chr"]).replace("chr", ""), int(d["pos"]))
    return pos


def main():
    os.makedirs(OUT, exist_ok=True)
    lift = ChainLift()
    pos37 = load_cpg_pos()
    print("CpGs with coordinates: %d / %d" % (len(pos37), len(CPGS)), flush=True)

    tb = TabixHTTP(URL, IDX)
    for region, cpgs in REGIONS.items():
        fp = os.path.join(OUT, "igan_%s.tsv" % region)
        if os.path.exists(fp) and os.path.getsize(fp) > 1000:
            print("[skip] %-20s cached" % region, flush=True)
            continue
        known = [(c, pos37[c]) for c in cpgs if c in pos37]
        if not known:
            print("[warn] %-20s no coordinates" % region, flush=True)
            continue
        chrom = known[0][1][0]
        coords38 = []
        for c, (ch, p37) in known:
            if ch != chrom:
                continue
            p38 = lift.map(ch, p37)
            if p38:
                coords38.append(p38)
        if not coords38:
            print("[warn] %-20s liftover failed" % region, flush=True)
            continue
        beg = max(min(coords38) - WINDOW, 1)
        end = max(coords38) + WINDOW
        try:
            lines = tb.query(chrom, beg, end, max_blocks=12000)
        except Exception as ex:                       # noqa: BLE001
            print("[FAIL] %-20s %s" % (region, ex), flush=True)
            continue
        out = []
        for f in lines:
            if len(f) < 12 or f[0].startswith("#") or f[0] == "chromosome":
                continue
            out.append({"rsid": f[11] if f[11] not in ("", "NA", ".") else "",
                        "chr": f[0], "pos38": f[1],
                        "ea": f[2], "oa": f[3],
                        "beta": f[4], "se": f[5], "eaf": f[6], "p": f[7]})
        with open(fp, "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=FIELDS)
            w.writeheader()
            w.writerows(out)
        n_id = sum(1 for r in out if r["rsid"])
        print("[ok]   %-20s chr%s:%d-%d  rows=%-7d with_rsid=%-7d cpgs=%d"
              % (region, chrom, beg, end, len(out), n_id, len(known)), flush=True)


if __name__ == "__main__":
    sys.exit(main())
