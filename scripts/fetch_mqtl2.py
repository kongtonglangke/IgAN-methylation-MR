"""Fetch GoDMC cis-mQTL for candidate regions (GRCh37).
Writes one JSON per gene into data/mqtl/ for resumability.
Avoids heredoc so no backslash mangling.
"""
import json, os, ssl, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

WIN = int(sys.argv[1]) if len(sys.argv) > 1 else 500000
OUTDIR = "data/mqtl"
os.makedirs(OUTDIR, exist_ok=True)

genes = []
with open("data/gene_coords_grch37.tsv") as f:
    next(f)
    for line in f:
        p = line.rstrip().split(chr(9))
        if len(p) < 4 or p[1] == "NA":
            continue
        genes.append((p[0], p[1], int(p[2]), int(p[3])))


def log(msg):
    print(msg, flush=True)


def fetch(g):
    name, chrom, start, end = g
    path = os.path.join(OUTDIR, name + ".json")
    if os.path.exists(path) and os.path.getsize(path) > 2:
        try:
            with open(path) as fh:
                rows = json.load(fh)
            return (name, len(rows), "cached")
        except Exception:
            pass
    lo = max(1, start - WIN)
    hi = end + WIN
    url = "http://api.godmc.org.uk/v0.1/assoc_meta/range/cpg/%s:%d-%d" % (chrom, lo, hi)
    last = ""
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=600, context=CTX) as r:
                chunks = []
                while True:
                    c = r.read(1 << 20)
                    if not c:
                        break
                    chunks.append(c)
            raw = b"".join(chunks).decode("utf-8", "ignore")
            d = json.loads(raw)
            rows = d.get("assoc_meta", [])
            with open(path, "w") as fh:
                json.dump(rows, fh)
            return (name, len(rows), "ok")
        except Exception as ex:
            last = type(ex).__name__ + ":" + str(ex)[:60]
            time.sleep(5 + attempt * 5)
    return (name, 0, "ERR " + last)


if __name__ == "__main__":
    log("genes=%d win=%d" % (len(genes), WIN))
    with ThreadPoolExecutor(max_workers=3) as ex:
        for name, n, st in ex.map(fetch, genes):
            log("%-12s rows=%-6d %s" % (name, n, st))
    total_files = len([x for x in os.listdir(OUTDIR) if x.endswith(".json")])
    log("DONE gene_files=%d" % total_files)
