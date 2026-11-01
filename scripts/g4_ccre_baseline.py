# -*- coding: utf-8 -*-
"""G4 补充: 随机窗口 cCRE 密度基线

在 chr8 与 chr22 上各取若干等宽窗口，计算 cCRE 密度分布，
用于判断两个目标区域是否显著偏离局部背景。
"""
import io
import json
import os
import random
import time
import urllib.request
import ssl

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
API = ("https://api.genome.ucsc.edu/getData/track?"
       "genome=hg38;track=encodeCcreCombined;chrom={chrom};start={s};end={e}")

CHROM_SIZE = {"chr8": 145138636, "chr22": 50818468}
WIN = {"chr8": 110000, "chr22": 40000}
TARGETS = {"chr8": (6920000, 7030000), "chr22": (30176000, 30216000)}
N_RANDOM = 12


def count_ccre(chrom, s, e, retries=3):
    for i in range(retries):
        try:
            req = urllib.request.Request(API.format(chrom=chrom, s=s, e=e),
                                         headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60, context=CTX) as r:
                d = json.loads(r.read().decode())
            return len(d.get("encodeCcreCombined", []))
        except Exception:
            time.sleep(2 * (i + 1))
    return None


def main():
    random.seed(20261101)
    rows = []
    for chrom in ("chr8", "chr22"):
        w = WIN[chrom]
        ts, te = TARGETS[chrom]
        n = count_ccre(chrom, ts, te)
        rows.append((chrom, "TARGET", ts, te, n, round(n / (w / 1000.0), 3)))
        print("%s TARGET %d-%d n=%s" % (chrom, ts, te, n), flush=True)
        got = 0
        tries = 0
        while got < N_RANDOM and tries < N_RANDOM * 6:
            tries += 1
            s = random.randint(0, CHROM_SIZE[chrom] - w - 1)
            e = s + w
            # 避免与目标窗口重叠
            if chrom == "chr8" and not (e < ts or s > te):
                continue
            if chrom == "chr22" and not (e < ts or s > te):
                continue
            n = count_ccre(chrom, s, e)
            if n is None:
                continue
            rows.append((chrom, "RANDOM", s, e, n, round(n / (w / 1000.0), 3)))
            got += 1
            print("  %s rand %d-%d n=%d" % (chrom, s, e, n), flush=True)

    outp = os.path.join(RES, "g4_ccre_density_baseline.tsv")
    with io.open(outp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("chrom\twindow_type\tstart\tend\tn_ccre\tccre_per_kb\n")
        for r in rows:
            fh.write("\t".join(str(v) for v in r) + "\n")
    print("写出", outp)

    # 汇总
    import statistics as st
    for chrom in ("chr8", "chr22"):
        tgt = [r for r in rows if r[0] == chrom and r[1] == "TARGET"][0]
        rnd = [r[5] for r in rows if r[0] == chrom and r[1] == "RANDOM"]
        if len(rnd) < 3:
            continue
        rnd_sorted = sorted(rnd)
        rank = sum(1 for v in rnd if v < tgt[5])
        print("%s: target %.3f/kb | random n=%d median %.3f (range %.3f-%.3f) | "
              "percentile %.0f%%" % (chrom, tgt[5], len(rnd), st.median(rnd),
                                     min(rnd), max(rnd), 100.0 * rank / len(rnd)))


if __name__ == "__main__":
    main()
