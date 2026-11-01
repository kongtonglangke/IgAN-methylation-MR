#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fetch_mqtldb_snp.py —— 按 SNP 查 mQTLdb(ARIES)，返回该 SNP 的全部 cis-mQTL CpG

用途：cg00705004 在 ARIES 里被 QC 剔除（0 行、且不在 CpG 位置表内），
无法直接复现。改为问一个区域性问题：
  「IgAN 风险单倍型在 ARIES 里是否仍是 8p23.1 defensin 簇内 CpG 的 mQTL？」
按 SNP 查询可一次拿到该 SNP 在该区域调用的所有 CpG，避免逐个 CpG 猜。

输出 results/aries_snp_<rsid>.tsv
"""
import os, re, sys, io, csv, time, subprocess, threading, queue

BASE = "http://www.mqtldb.org/cgi-bin/search.cgi"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAWD = os.path.join(ROOT, "data", "mqtldb", "raw_snp")
OUTD = os.path.join(ROOT, "results")
os.makedirs(RAWD, exist_ok=True)

SNPS = sys.argv[1:] or [
    "rs2738048", "rs4512398", "rs12716641", "rs10086568",   # DEFA1 8p23.1 单倍型
    "rs73164728", "rs12537",                                # HORMAD2 22q12
    "rs6677604",                                            # CFH 1q32
]
COLS = ["timepoint", "snp", "snp_chr", "snp_pos", "a1", "a2", "maf",
        "cpg", "cpg_chr", "cpg_pos", "beta", "tstat", "effect_size", "pval", "trans"]
LOCK = threading.Lock()


def fetch(rsid, attempts=5):
    raw = os.path.join(RAWD, "%s.html" % rsid)
    for a in range(1, attempts + 1):
        tmp = raw + ".tmp"
        cmd = ["curl", "-s", "--ssl-no-revoke", "--max-time", "300",
               "--speed-limit", "2000", "--speed-time", "25",
               "--form-string", "batch=%s" % rsid,
               "-F", "database=1", "-F", "timepoint=0", "-F", "distance=1000000",
               BASE, "-o", tmp]
        subprocess.run(cmd, capture_output=True)
        ok = False
        if os.path.exists(tmp):
            with io.open(tmp, "r", encoding="utf-8", errors="replace") as fh:
                ok = "</table>" in fh.read()          # 先关句柄再改名，否则 Windows 会报占用
            if ok:
                os.replace(tmp, raw)
                return True
            os.remove(tmp)
        with LOCK:
            print("  [retry] %s a=%d" % (rsid, a), flush=True)
        time.sleep(10 * a)
    return False


def parse(path):
    with io.open(path, "r", encoding="utf-8", errors="replace") as fh:
        html = fh.read()
    out = []
    for tr in re.findall(r"<tr>(.*?)</tr>", html, flags=re.S):
        tds = [re.sub(r"<[^>]+>", "", x).strip() for x in
               re.findall(r"<td[^>]*>(.*?)</td>", tr, flags=re.S)]
        if len(tds) >= 14:
            out.append(tds[:15] if len(tds) >= 15 else tds + [""] * (15 - len(tds)))
    return out


def worker(q):
    while True:
        try:
            rsid = q.get_nowait()
        except queue.Empty:
            return
        raw = os.path.join(RAWD, "%s.html" % rsid)
        if not (os.path.exists(raw) and os.path.getsize(raw) > 5000):
            if not fetch(rsid):
                print("  !! %s 失败" % rsid, flush=True)
                q.task_done()
                continue
        rows = [r for r in parse(raw) if r[1] == rsid]
        outp = os.path.join(OUTD, "aries_snp_%s.tsv" % rsid)
        with io.open(outp, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\t".join(COLS) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        print("  -> %-12s %4d CpG" % (rsid, len(rows)), flush=True)
        q.task_done()


def main():
    q = queue.Queue()
    for s in SNPS:
        q.put(s)
    ths = [threading.Thread(target=worker, args=(q,), daemon=True) for _ in range(4)]
    for t in ths:
        t.start()
    for t in ths:
        t.join()


if __name__ == "__main__":
    main()
