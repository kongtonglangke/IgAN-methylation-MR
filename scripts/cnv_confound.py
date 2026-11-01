"""DEFA1 CNV-confounding diagnostics.

Three questions, all answerable from public reference data plus the mQTL
summary statistics we already hold:

  Q1  Is the CpG (or its 450K probe footprint) inside a segmental
      duplication?  -> UCSC genomicSuperDups, fracMatch + length
  Q2  How many near-identical copies of the probe-length sequence exist
      inside the defensin cluster?  -> local k-mer near-match search
  Q3  Does the lead SNP behave like a copy-number proxy (many CpGs shifted
      together, high between-cohort heterogeneity)?  -> GoDMC region dump

Q2 is done locally so it does not depend on BLAT, which is bot-protected.
"""
import json
import os
import ssl
import sys
import time
import urllib.request
from collections import Counter

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

API = "https://api.genome.ucsc.edu/"
CACHE = "data/ref"
PROBE_LEN = 50           # Infinium probe is 50 nt


def api(path, timeout=120, tries=3):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(API + path, headers={"User-Agent": "research/1.0"})
            with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
                return json.loads(r.read().decode("utf8", "ignore"))
        except Exception as ex:                        # noqa: BLE001
            last = ex
            time.sleep(2 + 2 * i)
    raise last


def get_seq(chrom, start, end, genome="hg19", step=20000):
    """Cached reference sequence, small chunks with resume so a dropped
    connection does not lose the whole download."""
    os.makedirs(CACHE, exist_ok=True)
    fp = os.path.join(CACHE, "%s_%s_%d_%d.txt" % (genome, chrom, start, end))
    if os.path.exists(fp) and len(open(fp).read()) == end - start:
        return open(fp).read()
    pos = start
    if os.path.exists(fp):
        pos = start + len(open(fp).read())
    fh = open(fp, "a")
    while pos < end:
        n = min(step, end - pos)
        for attempt in range(5):
            try:
                d = api("getData/sequence?genome=%s;chrom=%s;start=%d;end=%d"
                        % (genome, chrom, pos, pos + n), timeout=60, tries=2)
                fh.write(d["dna"])
                fh.flush()
                pos += n
                break
            except Exception as ex:                    # noqa: BLE001
                print("    retry seq at %d (%s)" % (pos, type(ex).__name__), flush=True)
                time.sleep(3 + 3 * attempt)
        else:
            raise RuntimeError("sequence fetch stalled at %d" % pos)
        time.sleep(0.25)
    fh.close()
    return open(fp).read()


RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(RC)[::-1]


def kmers(seq, k):
    seq = seq.upper()
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]


def near_match_count(region, query, max_mm):
    """Count occurrences of `query` in `region` allowing <= max_mm mismatches.

    Only exact-length windows are compared, which is what a 50-mer probe does.
    """
    q = query.upper()
    k = len(q)
    hits = []
    # seed on the first 12 nt, then verify — much faster than brute force
    seed = q[:12]
    start = 0
    while True:
        i = region.find(seed, start)
        if i < 0:
            break
        start = i + 1
        cand = region[i:i + k]
        if len(cand) < k:
            continue
        mm = sum(1 for a, b in zip(cand, q) if a != b)
        if mm <= max_mm:
            hits.append((i, mm))
    return hits


def load_segdups(chrom, start, end):
    d = api("getData/track?genome=hg19;track=genomicSuperDups;chrom=%s;start=%d;end=%d"
            % (chrom, start, end))
    return d.get("genomicSuperDups", []) or []


def gene_annotation(chrom, start, end):
    d = api("getData/track?genome=hg19;track=ncbiRefSeqCurated;chrom=%s;start=%d;end=%d"
            % (chrom, start, end))
    return d.get("ncbiRefSeqCurated", []) or []


def main():
    cpgs = {
        # cpg            gene        chrom  pos37      role
        "cg00705004":  ("DEFA1",     "chr8",  6827814, "coloc H4=0.976"),
        "cg07297068":  ("DEFA1",     "chr8",  6691691, "same region"),
        "cg19304352":  ("DEFA1",     "chr8",  6796176, "same region"),
        "cg05969298":  ("DEFA1",     "chr8",  6924737, "same region"),
        "cg00539347":  ("HORMAD2",   "chr22", None,    "coloc H4=0.982 (control)"),
        "cg13682187":  ("CFH",       "chr1",  None,    "coloc H3 (control)"),
    }
    # fill the control coordinates from the GoDMC info cache
    info = json.load(open("data/coloc/cpg_info.json"))
    for c, v in info.items():
        if c in cpgs and cpgs[c][2] is None:
            x = v[0] if isinstance(v, list) else v
            cpgs[c] = (cpgs[c][0], "chr" + str(x["chr"]), int(x["pos"]), cpgs[c][3])

    out = {"segdup": {}, "kmer": {}, "annotation": {}}

    print("=" * 78)
    print("Q1  segmental duplications overlapping each CpG (hg19, 19kb+ identity tiles)")
    print("=" * 78)
    seg = load_segdups("chr8", 6500000, 7200000)
    print("records in chr8:6.5-7.2Mb: %d  (keeping len>=2kb and fracMatch>=0.90)" % len(seg))
    blocks = []
    for r in seg:
        if r["chromEnd"] - r["chromStart"] < 2000 or r["fracMatch"] < 0.90:
            continue
        blocks.append((r["chromStart"], r["chromEnd"], r["fracMatch"],
                       r["otherChrom"], r["otherStart"], r["otherEnd"]))
    # collapse the reciprocal pairs into unique intervals
    uniq = sorted({(b[0], b[1], round(b[2], 4)) for b in blocks} |
                  {(b[4], b[5], round(b[2], 4)) for b in blocks})
    print("%-24s %8s %10s" % ("interval (hg19)", "len", "fracMatch"))
    for a, b, f in uniq:
        print("%-24s %8d %10.4f" % ("%d-%d" % (a, b), b - a, f))
    out["segdup"]["chr8_6.5_7.2"] = uniq

    print("\nCpG position vs those intervals:")
    for c, (gene, chrom, pos, role) in cpgs.items():
        if chrom != "chr8" or pos is None:
            continue
        inside = [(a, b, f) for a, b, f in uniq if a <= pos < b]
        near = min((min(abs(pos - a), abs(pos - b)) for a, b, f in uniq), default=None)
        print("  %-12s %-9s %8d  inside=%s  nearest_edge=%s bp" %
              (c, gene, pos, ("YES frac=%.4f" % inside[0][2]) if inside else "no", near))
        out["segdup"][c] = dict(gene=gene, pos=pos, inside=bool(inside),
                                frac=inside[0][2] if inside else None, nearest=near)

    # ---------------------------------------------------------------- Q2
    print("\n" + "=" * 78)
    print("Q2  near-identical copies of a %d nt probe footprint in chr8:6.78-6.95Mb" % PROBE_LEN)
    print("=" * 78)
    region = get_seq("chr8", 6780000, 6950000)
    print("region downloaded: %d bp" % len(region))
    for c, (gene, chrom, pos, role) in cpgs.items():
        if chrom != "chr8" or pos is None:
            continue
        rel = pos - 6780000
        # probe sits on either side of the CpG; test both orientations
        for lab, s, e in (("3'", rel, rel + PROBE_LEN), ("5'", rel - PROBE_LEN, rel)):
            if s < 0 or e > len(region):
                continue
            q = region[s:e]
            for strand, qq in (("+", q), ("-", revcomp(q))):
                hits = near_match_count(region, qq, max_mm=2)
                # drop the self hit
                hits = [h for h in hits if abs(h[0] - s) > 3]
                if hits:
                    print("  %-12s %-9s probe%s/%s  extra near-copies(<=2mm): %d  at %s"
                          % (c, gene, lab, strand, len(hits),
                             [("%d(mm=%d)" % (6780000 + h[0], h[1])) for h in hits[:4]]))
                else:
                    print("  %-12s %-9s probe%s/%s  extra near-copies(<=2mm): 0"
                          % (c, gene, lab, strand))
        out["kmer"][c] = True

    print("\n" + "=" * 78)
    print("Q3  defensin cluster gene models (hg19 6.7-7.0Mb)")
    print("=" * 78)
    g = gene_annotation("chr8", 6700000, 7000000)
    seen = set()
    for r in sorted(g, key=lambda x: x["txStart"]):
        if r["name2"] in seen:
            continue
        seen.add(r["name2"])
        print("  %-10s %8d-%8d %s" % (r["name2"], r["txStart"], r["txEnd"], r["strand"]))

    json.dump(out, open("results/cnv_confound_ref.json", "w"), indent=1)


if __name__ == "__main__":
    main()
