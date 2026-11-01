#!/usr/bin/env python
"""Fetch FinnGen R12 IgA nephropathy summary statistics for the coloc regions.

Why FinnGen and not Sakaue: GCST90018866 is the only IgAN GWAS with public
summary statistics, but its signal is dominated by the MHC — in the CFH,
HORMAD2, ITGAM, MADCAM1, ST6GAL1 and TNFRSF13B regions no SNP reaches
p<5e-8, which makes colocalisation impossible.  FinnGen R12 N14_IGA_NEPHROPATHY
is a registry/ biopsy based endpoint with full summary statistics and a tabix
index.

FinnGen columns:
    #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta
    af_alt af_alt_cases af_alt_controls
Coordinates are GRCh38, so the GRCh37 CpG positions are lifted once.

Output: data/coloc/finngen_<region>.tsv
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from liftover_chain import ChainLift          # noqa: E402
from tabix_http import TabixHTTP             # noqa: E402

BASE = ("https://storage.googleapis.com/finngen-public-data-r12/"
        "summary_stats/release/finngen_R12_N14_IGA_NEPHROPATHY.gz")
URL_GZ = BASE
URL_TBI = BASE + ".tbi"
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


def load_coords():
    """CpG positions in GRCh37 from the cached GoDMC info payload."""
    info = json.load(open(os.path.join(OUT, "cpg_info.json")))
    out = {}
    for cpg in {c for v in REGIONS.values() for c in v}:
        if cpg not in info:
            continue
        d = info[cpg]
        d = d[0] if isinstance(d, list) else d
        if d.get("chr") and d.get("pos"):
            out[cpg] = (str(d["chr"]).replace("chr", ""), int(d["pos"]))
    return out


def region_window(region):
    """Reuse the exact window that was queried in the Sakaue file.

    Identical windows on both sides make the SNP intersection cleaner and
    avoids depending on the CpG coordinate cache (some CpGs were reached from
    the IgAN side and never got an entry there).
    """
    fp = os.path.join(OUT, "igan_%s.tsv" % region)
    if not os.path.exists(fp):
        return None
    pos = []
    chrom = None
    for r in csv.DictReader(open(fp), delimiter="\t"):
        v = r.get("pos38", "")
        if v.isdigit():
            pos.append(int(v))
        if chrom is None and str(r.get("chr", "")).isdigit():
            chrom = str(r["chr"])
    if not pos or chrom is None:
        return None
    return (min(pos), max(pos), chrom)


def cpG_window(region, coords, lift):
    chrom = None
    c38 = []
    for c in REGIONS[region]:
        if c not in coords:
            continue
        ch, p37 = coords[c]
        if chrom is None:
            chrom = ch
        if ch != chrom:
            continue
        p = lift.map(ch, p37)
        if p:
            c38.append(p)
    if not c38:
        return None
    return (max(min(c38) - WINDOW, 1), max(c38) + WINDOW, chrom)


def main():
    os.makedirs(OUT, exist_ok=True)
    lift = ChainLift()
    coords = load_coords()
    print("CpG coordinates available: %d" % len(coords), flush=True)

    tb = TabixHTTP(URL_GZ, URL_TBI, cache_dir="data/tbi_cache_finngen")
    for region, cpgs in REGIONS.items():
        fp = os.path.join(OUT, "finngen_%s.tsv" % region)
        if os.path.exists(fp) and os.path.getsize(fp) > 1000:
            print("[skip] %-20s cached" % region, flush=True)
            continue

        win = region_window(region)
        if win:
            beg, end, chrom = win
            print("       %-20s window from Sakaue file chr%s:%d-%d"
                  % (region, chrom, beg, end), flush=True)
        else:
            w = cpG_window(region, coords, lift)
            if not w:
                print("[warn] %-20s no window available" % region, flush=True)
                continue
            beg, end, chrom = w
        try:
            lines = tb.query(chrom, beg, end, max_blocks=12000)
        except Exception as ex:                       # noqa: BLE001
            print("[FAIL] %-20s %s" % (region, ex), flush=True)
            continue
        out = []
        for f in lines:
            if len(f) < 11 or f[0].startswith("#"):
                continue
            rsid = f[4].split(",")[0] if f[4] not in ("", "NA", ".") else ""
            out.append({"rsid": rsid, "chr": f[0], "pos38": f[1],
                        "ea": f[3], "oa": f[2],
                        "beta": f[8], "se": f[9], "eaf": f[10], "p": f[6]})
        with open(fp, "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=FIELDS)
            w.writeheader()
            w.writerows(out)
        ps = sorted(float(r["p"]) for r in out if r["p"] not in ("", "NA"))
        print("[ok]   %-20s chr%s:%d-%d rows=%-6d with_rsid=%-6d min_p=%.2e p<5e-8=%d"
              % (region, chrom, beg, end, len(out),
                 sum(1 for r in out if r["rsid"]),
                 ps[0] if ps else 0, sum(1 for p in ps if p < 5e-8)), flush=True)


if __name__ == "__main__":
    sys.exit(main())
