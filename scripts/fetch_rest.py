"""Sequential (throttled) fetch of remaining GoDMC regions.
GoDMC rate-limits aggressively: one region at a time, with a pause between calls.
Skips genes already cached in data/mqtl/.
"""
import json, os, ssl, time, urllib.request

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

WIN = 500000
DELAY = 15          # seconds between successful calls
OUTDIR = "data/mqtl"
os.makedirs(OUTDIR, exist_ok=True)

genes = []
with open("data/gene_coords_grch37.tsv") as f:
    next(f)
    for line in f:
        p = line.rstrip().split(chr(9))
        if len(p) >= 4 and p[1] != "NA":
            genes.append((p[0], p[1], int(p[2]), int(p[3])))

todo = []
for g in genes:
    path = os.path.join(OUTDIR, g[0] + ".json")
    if os.path.exists(path) and os.path.getsize(path) > 2:
        continue
    todo.append(g)

print("todo=%d" % len(todo), flush=True)
for i, (name, chrom, start, end) in enumerate(todo):
    lo = max(1, start - WIN)
    hi = end + WIN
    url = "http://api.godmc.org.uk/v0.1/assoc_meta/range/cpg/%s:%d-%d" % (chrom, lo, hi)
    ok = False
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=600, context=CTX) as r:
                buf = []
                while True:
                    c = r.read(1 << 20)
                    if not c:
                        break
                    buf.append(c)
            rows = json.loads(b"".join(buf).decode("utf-8", "ignore")).get("assoc_meta", [])
            with open(os.path.join(OUTDIR, name + ".json"), "w") as fh:
                json.dump(rows, fh)
            print("[%d/%d] %-12s rows=%-6d OK" % (i + 1, len(todo), name, len(rows)), flush=True)
            ok = True
            break
        except Exception as ex:
            code = getattr(ex, "code", None)
            print("[%d/%d] %-12s attempt %d ERR %s" % (i + 1, len(todo), name, attempt + 1, str(ex)[:60]), flush=True)
            if code == 403:
                time.sleep(60 * (attempt + 1))
            else:
                time.sleep(10 * (attempt + 1))
    if not ok:
        print("[%d/%d] %-12s FAILED" % (i + 1, len(todo), name), flush=True)
    time.sleep(DELAY)
print("DONE", flush=True)
