#!/usr/bin/env python
"""Fetch full cis-mQTL records for the CpGs that colocalise with IgAN risk loci.

GoDMC's public database only stores the candidate set (the union of SNP-CpG
pairs that reached p<1e-5 in at least one cohort), so this is the densest
mQTL data that exists publicly.  The POST /v0.1/query endpoint with no pval
field applies no additional filter.

Output: data/coloc/mqtl_<cpg>.json  and data/coloc/cpg_info.json
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

# one entry per independent IgAN risk region; a few CpGs each so that the
# strongest methylation signal per region is represented
TARGETS = {
    "CFH_CFHR_1q32":    ["cg13682187", "cg11348583", "cg13441112", "cg07209298"],
    "DEFA1_8p23":       ["cg00705004", "cg07297068", "cg19304352", "cg05969298"],
    "ITGAM_ITGAX_16p11": ["cg03100639", "cg10451425", "cg02846316", "cg04018474"],
    "HORMAD2_22q12":    ["cg27665648", "cg00539347"],
    "MADCAM1_19p13":    ["cg23386659", "cg27535410"],
    "ST6GAL1_3q27":     ["cg17448192", "cg25599673"],
    "TNFRSF13B_17p11":  ["cg02671941"],
}


def post(body, timeout=600, tries=3):
    data = json.dumps(body).encode()
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(
                API + "/query", data=data,
                headers={"Content-Type": "application/json",
                         "Accept": "application/json",
                         "User-Agent": "igan-epi-mr/1.0"})
            with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as ex:                       # noqa: BLE001
            last = ex
            time.sleep(3 + 4 * i)
    raise RuntimeError(str(last))


def get_json(path, timeout=120, tries=3):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(API + path,
                                         headers={"Accept": "application/json",
                                                  "User-Agent": "igan-epi-mr/1.0"})
            with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as ex:                       # noqa: BLE001
            last = ex
            time.sleep(2 + 3 * i)
    return {"error": str(last)}


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
    info_path = os.path.join(OUT, "cpg_info.json")
    info = json.load(open(info_path)) if os.path.exists(info_path) else {}

    for region, cpgs in TARGETS.items():
        for cpg in cpgs:
            fp = os.path.join(OUT, "mqtl_%s.json" % cpg)
            if os.path.exists(fp) and os.path.getsize(fp) > 100:
                n = len(json.load(open(fp)))
                print("[skip] %-12s %-16s cached n=%d" % (cpg, region, n), flush=True)
                continue
            t0 = time.time()
            try:
                arr = unwrap(post({"cpgs": [cpg], "cistrans": "cis"}))
            except Exception as ex:                   # noqa: BLE001
                print("[FAIL] %-12s %s" % (cpg, ex), flush=True)
                continue
            json.dump(arr, open(fp, "w"))
            ps = [x.get("pval") for x in arr if x.get("pval") is not None]
            n8 = sum(1 for p in ps if p < 5e-8)
            print("[ok]   %-12s %-16s n=%-5d p<5e-8=%-4d maxp=%.1e  %.0fs"
                  % (cpg, region, len(arr), n8, max(ps) if ps else 0,
                     time.time() - t0), flush=True)

            if cpg not in info:
                try:
                    info[cpg] = get_json("/info/cpg/%s" % cpg)
                    json.dump(info, open(info_path, "w"), indent=1)
                except Exception:                     # noqa: BLE001
                    pass
            time.sleep(1.5)

    print("\n=== summary ===")
    tot = 0
    for region, cpgs in TARGETS.items():
        ns = []
        for cpg in cpgs:
            fp = os.path.join(OUT, "mqtl_%s.json" % cpg)
            if os.path.exists(fp):
                ns.append((cpg, len(json.load(open(fp)))))
                tot += ns[-1][1]
        print("%-20s %s" % (region, ", ".join("%s:%d" % x for x in ns)))
    print("total mQTL rows: %d" % tot)


if __name__ == "__main__":
    sys.exit(main())
