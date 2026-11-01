#!/usr/bin/env python
"""Retry the CpGs that the GoDMC API refused with 403 during the bulk pull.

The API starts returning 403 after a burst of requests, so this version waits
between calls and backs off exponentially when refused.  It only touches CpGs
whose output file is still missing.

Output: data/coloc/mqtl_<cpg>.json
"""
import json
import os
import ssl
import sys
import time
import urllib.request

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

API = "http://api.godmc.org.uk/v0.1"
OUT = "data/coloc"

PENDING = {
    "cg10451425": "ITGAM_ITGAX_16p11",
    "cg02846316": "ITGAM_ITGAX_16p11",
    "cg04018474": "ITGAM_ITGAX_16p11",
    "cg27665648": "HORMAD2_22q12",
    "cg00539347": "HORMAD2_22q12",
    "cg23386659": "MADCAM1_19p13",
    "cg27535410": "MADCAM1_19p13",
    "cg17448192": "ST6GAL1_3q27",
    "cg25599673": "ST6GAL1_3q27",
    "cg02671941": "TNFRSF13B_17p11",
}

GAP = 40          # seconds between successful calls
MAX_BACKOFF = 600


def post(body, timeout=600):
    req = urllib.request.Request(
        API + "/query", data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json",
                 "Accept": "application/json",
                 "User-Agent": "igan-epi-mr/1.0"})
    with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
        return json.loads(r.read().decode("utf-8"))


def unwrap(js):
    if isinstance(js, list):
        return js
    if isinstance(js, dict):
        for k in ("assoc_meta", "data", "results"):
            if k in js:
                v = js[k]
                return v if isinstance(v, list) else list(v.values())
        return list(js.values())
    return []


def main():
    os.makedirs(OUT, exist_ok=True)
    done, failed = [], []
    for cpg, region in PENDING.items():
        fp = os.path.join(OUT, "mqtl_%s.json" % cpg)
        if os.path.exists(fp) and os.path.getsize(fp) > 100:
            print("[skip] %-12s cached" % cpg, flush=True)
            continue

        backoff = 90
        got = None
        for attempt in range(1, 6):
            t0 = time.time()
            try:
                arr = unwrap(post({"cpgs": [cpg], "cistrans": "cis"}))
                if not arr:
                    print("[warn] %-12s empty payload (attempt %d)" % (cpg, attempt),
                          flush=True)
                got = arr
                break
            except Exception as ex:                   # noqa: BLE001
                msg = str(ex)
                print("[%s] %-12s attempt %d failed: %s" %
                      ("403" if "403" in msg else "ERR", cpg, attempt, msg),
                      flush=True)
                if attempt == 5:
                    break
                wait = min(backoff, MAX_BACKOFF)
                print("        sleeping %ds before retry" % wait, flush=True)
                time.sleep(wait)
                backoff = int(backoff * 1.7)

        if got is None:
            failed.append(cpg)
            continue
        json.dump(got, open(fp, "w"))
        ps = [x.get("pval") for x in got if x.get("pval") is not None]
        n8 = sum(1 for p in ps if p < 5e-8)
        print("[ok]   %-12s %-18s n=%-5d p<5e-8=%-4d  %.0fs"
              % (cpg, region, len(got), n8, time.time() - t0), flush=True)
        done.append(cpg)
        time.sleep(GAP)

    print("\nrecovered: %d  failed: %d %s" % (len(done), len(failed), failed))
    have = [f for f in os.listdir(OUT) if f.startswith("mqtl_")]
    print("mqtl files now: %d" % len(have))


if __name__ == "__main__":
    sys.exit(main())
