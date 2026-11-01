"""Fetch GoDMC mQTL in small chunks to avoid the server's large-query 403 limit.
Each gene window is split into CHUNK-sized pieces; results are de-duplicated on
(cpg, snp) and written to data/mqtl/<gene>.json in the same schema as the
whole-window fetch, so downstream code needs no change.
"""
import json, os, ssl, time, urllib.request

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

WIN = 500000
CHUNK = 100000
DELAY = 4
OUTDIR = "data/mqtl"
os.makedirs(OUTDIR, exist_ok=True)


def fetch_range(chrom, lo, hi):
    url = "http://api.godmc.org.uk/v0.1/assoc_meta/range/cpg/%s:%d-%d" % (chrom, lo, hi)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=300, context=CTX) as r:
                buf = []
                while True:
                    c = r.read(1 << 20)
                    if not c:
                        break
                    buf.append(c)
            return json.loads(b"".join(buf).decode("utf-8", "ignore")).get("assoc_meta", [])
        except Exception as ex:
            code = getattr(ex, "code", None)
            wait = 60 * (attempt + 1) if code == 403 else 10 * (attempt + 1)
            print("     chunk %s:%d-%d attempt %d %s (wait %ds)" % (chrom, lo, hi, attempt + 1, str(ex)[:40], wait), flush=True)
            time.sleep(wait)
    return None


def main():
    genes = []
    with open("data/gene_coords_grch37.tsv") as f:
        next(f)
        for line in f:
            p = line.rstrip().split(chr(9))
            if len(p) >= 4 and p[1] != "NA":
                genes.append((p[0], p[1], int(p[2]), int(p[3])))

    todo = [g for g in genes if not (os.path.exists(os.path.join(OUTDIR, g[0] + ".json"))
                                     and os.path.getsize(os.path.join(OUTDIR, g[0] + ".json")) > 2)]
    print("todo=%d of %d" % (len(todo), len(genes)), flush=True)

    for gi, (name, chrom, start, end) in enumerate(todo):
        lo = max(1, start - WIN)
        hi = end + WIN
        merged = {}
        ok_all = True
        cur = lo
        nchunk = 0
        while cur <= hi:
            cend = min(hi, cur + CHUNK - 1)
            rows = fetch_range(chrom, cur, cend)
            nchunk += 1
            if rows is None:
                ok_all = False
                print("[%d/%d] %-10s chunk %d FAILED, abort gene" % (gi + 1, len(todo), name, nchunk), flush=True)
                break
            for r in rows:
                merged[(r.get("cpg"), r.get("snp"))] = r
            cur = cend + 1
            time.sleep(DELAY)
        if ok_all:
            out = list(merged.values())
            with open(os.path.join(OUTDIR, name + ".json"), "w") as fh:
                json.dump(out, fh)
            print("[%d/%d] %-10s chunks=%-3d rows=%d OK" % (gi + 1, len(todo), name, nchunk, len(out)), flush=True)
        else:
            print("[%d/%d] %-10s SKIPPED (server limit)" % (gi + 1, len(todo), name), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
