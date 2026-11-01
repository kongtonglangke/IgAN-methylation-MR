#!/usr/bin/env python
"""Fetch every GWAS Catalog association for IgA nephropathy (MONDO_0005342)
and test whether the known IgAN risk loci overlap the 27 metabolite-target
genes that the mQTL instrument set was built on.

This is a feasibility gate for option A ("are IgAN risk loci mQTLs of the
target genes?").  If the two sets do not intersect at all, option A cannot
be answered by annotation -- the question itself would be empty.

Outputs
-------
data/gwas/igan_catalog_assoc.tsv   one row per association, deduplicated by locus
data/gwas/igan_loci_overlap.tsv    overlap table vs the 27 target genes
"""
import csv
import json
import os
import ssl
import sys
import time
import urllib.request

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

BASE = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
EFO = "MONDO_0005342"          # IgA glomerulonephritis
WINDOW = 500_000               # same +-500 kb window used for the mQTL pull


def get(url, timeout=180, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(
                url, headers={"Accept": "application/json",
                              "User-Agent": "igan-epi-mr/1.0"})
            with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as ex:                       # noqa: BLE001
            last = ex
            time.sleep(2 + 3 * i)
    raise RuntimeError("failed %s: %s" % (url, last))


def fetch_all(size=100):
    out, page = [], 0
    while True:
        url = "%s?efo_id=%s&size=%d&page=%d" % (BASE, EFO, size, page)
        js = get(url)
        recs = js.get("_embedded", {}).get("associations", [])
        total = js.get("page", {}).get("totalElements")
        out.extend(recs)
        print("  page %d: +%d (have %d / %s)" % (page, len(recs), len(out), total),
              flush=True)
        if not recs or total is None or len(out) >= int(total):
            break
        page += 1
        time.sleep(0.5)
    return out


def load_target_genes():
    """Gene coordinates (GRCh37) for the 27 metabolite-target genes."""
    path = "data/gene_coords_grch37.tsv"
    genes = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            chrom = r.get("chr") or r.get("chrom") or r.get("chromosome")
            sym = r.get("gene") or r.get("symbol")
            if not (sym and chrom and r.get("start")):
                continue
            genes[sym] = (str(chrom).replace("chr", ""),
                          int(float(r["start"])), int(float(r["end"])))
    return genes


def main():
    os.makedirs("data/gwas", exist_ok=True)
    print("fetching IgAN associations for %s ..." % EFO, flush=True)
    recs = fetch_all()
    print("total associations fetched: %d" % len(recs), flush=True)

    rows = []
    for a in recs:
        rsids = [s.get("rs_id") for s in (a.get("snp_allele") or []) if s.get("rs_id")]
        eff = [s.get("effect_allele") for s in (a.get("snp_allele") or [])]
        locs = a.get("locations") or []
        loc = locs[0] if locs else ""
        chrom = pos = ""
        if ":" in loc:
            chrom, pos = loc.split(":")[:2]
        rows.append({
            "rsid": rsids[0] if rsids else "",
            "effect_allele": eff[0] if eff else "",
            "loc38": loc,
            "chrom": chrom,
            "pos38": pos,
            "p_value": a.get("p_value"),
            "or_value": a.get("or_value"),
            "beta": a.get("beta"),
            "risk_freq": a.get("risk_frequency"),
            "mapped_genes": ",".join(a.get("mapped_genes") or []),
            "reported_trait": ",".join(a.get("reported_trait") or []),
            "accession": a.get("accession_id"),
            "pmid": a.get("pubmed_id"),
            "first_author": a.get("first_author"),
        })

    # dedupe by (loc, rsid)
    seen, uniq = set(), []
    for r in rows:
        k = (r["loc38"], r["rsid"])
        if k in seen:
            continue
        seen.add(k)
        uniq.append(r)
    uniq.sort(key=lambda r: (r["chrom"] or "99", int(r["pos38"] or 0)))
    print("unique locus-rsid rows: %d" % len(uniq), flush=True)

    with open("data/gwas/igan_catalog_assoc.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(uniq[0].keys()))
        w.writeheader()
        w.writerows(uniq)

    # ---- overlap with the 27 target genes -------------------------------
    genes = load_target_genes()
    print("\ntarget genes loaded: %d" % len(genes), flush=True)

    ov, gene_hits, span_hits = [], set(), 0
    for r in uniq:
        if not r["pos38"]:
            continue
        chrom, pos = str(r["chrom"]), int(r["pos38"])
        named = {g for g in (r["mapped_genes"] or "").split(",") if g}
        hit_named = sorted(named & set(genes))
        near = []
        for g, (gc, gs, ge) in genes.items():
            if gc != chrom:
                continue
            # coordinates are GRCh37; allow a generous 520 kb slack for the
            # assembly shift before flagging a candidate
            if gs - WINDOW - 20_000 <= pos <= ge + WINDOW + 20_000:
                d = 0 if gs <= pos <= ge else min(abs(pos - gs), abs(pos - ge))
                near.append((g, d, gs <= pos <= ge))
        if hit_named or near:
            span_hits += 1 if near else 0
            gene_hits.update(hit_named)
            for g, d, inside in near:
                ov.append({
                    "rsid": r["rsid"], "loc38": r["loc38"], "p_value": r["p_value"],
                    "or_value": r["or_value"], "mapped_genes": r["mapped_genes"],
                    "target_gene": g, "distance_bp": d,
                    "inside_gene": "yes" if inside else "no",
                    "named_in_catalog": "yes" if g in named else "no",
                    "accession": r["accession"], "pmid": r["pmid"],
                })

    ov.sort(key=lambda r: (r["target_gene"], r["distance_bp"]))
    with open("data/gwas/igan_loci_overlap.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["rsid", "loc38", "p_value", "or_value",
                                       "mapped_genes", "target_gene", "distance_bp",
                                       "inside_gene", "named_in_catalog",
                                       "accession", "pmid"])
        w.writeheader()
        w.writerows(ov)

    print("\n=== FEASIBILITY GATE ===")
    print("IgAN loci (unique)            : %d" % len(uniq))
    print("loci within +-500kb of a target gene : %d" % span_hits)
    print("target genes named by GWAS Catalog   : %d %s"
          % (len(gene_hits), sorted(gene_hits)))
    print("genes implicated by span (count)     :")
    from collections import Counter
    c = Counter(r["target_gene"] for r in ov)
    for g, n in c.most_common():
        best = min((r for r in ov if r["target_gene"] == g),
                   key=lambda r: r["distance_bp"])
        print("   %-10s n=%-3d nearest=%-10s %8s bp  p=%-10s OR=%s"
              % (g, n, best["rsid"], best["distance_bp"],
                 best["p_value"], best["or_value"]))


if __name__ == "__main__":
    sys.exit(main())
