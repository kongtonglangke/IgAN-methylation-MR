"""Scan Orru 2020 immune-cell trait metadata on EBI FTP.
Parses YAML line-by-line (no regex) to avoid shell-escaping issues.
Outputs data/orru_traits.tsv with gwas_id and trait name.
"""
import os, ssl, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"


def bucket(n):
    if 90001001 <= n <= 90002000:
        return "GCST90001001-GCST90002000"
    if 90002001 <= n <= 90003000:
        return "GCST90002001-GCST90003000"
    return None


def parse_trait(text):
    lines = text.split(chr(10))
    for i, ln in enumerate(lines):
        if ln.strip().startswith("trait_description:"):
            if i + 1 < len(lines):
                nxt = lines[i + 1].strip()
                if nxt.startswith("-"):
                    return nxt[1:].strip()
                return nxt
    return ""


def get(n):
    g = "GCST%d" % n
    b = bucket(n)
    if not b:
        return (g, "")
    url = BASE + b + "/" + g + "/" + g + "_buildGRCh37.tsv.gz-meta.yaml"
    for a in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30, context=CTX) as r:
                t = r.read().decode("utf-8", "ignore")
            tr = parse_trait(t)
            return (g, tr)
        except Exception:
            if a == 2:
                return (g, "")
            time.sleep(1)
    return (g, "")


if __name__ == "__main__":
    res = []
    with ThreadPoolExecutor(max_workers=8) as ex:
        for g, t in ex.map(get, range(90001391, 90002122)):
            if t:
                res.append((g, t))
    print("resolved=%d" % len(res), flush=True)
    os.makedirs("data", exist_ok=True)
    with open("data/orru_traits.tsv", "w") as f:
        f.write("gwas_id" + chr(9) + "trait" + chr(10))
        for g, t in res:
            f.write(g + chr(9) + t + chr(10))
    bc = [(g, t) for g, t in res if "B cell" in t and "T cell" not in t]
    print("B_cell_traits=%d" % len(bc), flush=True)
    for g, t in bc:
        print(g + chr(9) + t, flush=True)
