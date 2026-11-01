#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fetch_mqtldb.py —— 从 mQTLdb (ARIES, Gaunt 2016 Genome Biol) 取指定 CpG 的 mQTL 结果

设计要点
--------
1. mQTLdb 的结果表是 MatrixEQTL 输出，**只保留 FDR<0.05 的配对**，因此：
   - 可做「方向 + 效应量一致性」的独立队列复现
   - 不能做 coloc（缺未显著 SNP，缺 varbeta）
2. 单 CpG 抓取：服务端很慢（~20 KB/s），批量 19 个 CpG 一次请求会截断，故逐个取 + 重试。
3. 输出 results/aries_mqtl_<cpg>.tsv 与合并表 results/aries_mqtl_all.tsv
"""
import os, re, sys, time, subprocess, io

BASE = "http://www.mqtldb.org/cgi-bin/search.cgi"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAWD = os.path.join(ROOT, "data", "mqtldb", "raw")
OUTD = os.path.join(ROOT, "results")
os.makedirs(RAWD, exist_ok=True)
os.makedirs(OUTD, exist_ok=True)

CPGS = [
    "cg07209298", "cg11348583", "cg13441112", "cg13682187",          # CFH/CFHR 1q32
    "cg00705004", "cg05969298", "cg07297068", "cg19304352",          # DEFA1 8p23
    "cg00539347", "cg27665648",                                       # HORMAD2 22q12
    "cg02846316", "cg03100639", "cg10451425", "cg04018474",           # ITGAM/ITGAX 16p11
    "cg27535410", "cg23386659",                                       # MADCAM1 19p13
    "cg25599673", "cg17448192",                                       # ST6GAL1 3q27
    "cg02671941",                                                     # TNFRSF13B 17p11
]
COLS = ["timepoint", "snp", "snp_chr", "snp_pos", "a1", "a2", "maf",
        "cpg", "cpg_chr", "cpg_pos", "beta", "tstat", "effect_size", "pval", "trans"]


def curl_post(cpg, out_path, max_time=600):
    """POST 一条查询；-s 静默，失败返回 False"""
    cmd = ["curl", "-s", "--ssl-no-revoke", "-m", str(max_time),
           "--form-string", "batch=%s" % cpg,
           "-F", "database=1", "-F", "timepoint=0", "-F", "distance=1000000",
           BASE, "-o", out_path]
    r = subprocess.run(cmd, capture_output=True)
    return r.returncode == 0


def parse_html(path):
    """从结果页抽 <td> 行；返回 list[list[str]]"""
    with io.open(path, "r", encoding="utf-8", errors="replace") as fh:
        html = fh.read()
    if "</table>" not in html:          # 截断
        return None, 0
    rows = []
    for tr in re.findall(r"<tr>(.*?)</tr>", html, flags=re.S):
        tds = [re.sub(r"<[^>]+>", "", x).strip() for x in
               re.findall(r"<td[^>]*>(.*?)</td>", tr, flags=re.S)]
        if len(tds) >= 14:
            rows.append(tds[:15] if len(tds) >= 15 else tds + [""] * (15 - len(tds)))
    return rows, len(rows)


def main():
    global CPGS
    if len(sys.argv) > 1:                     # 可只抓指定 CpG： python fetch_mqtldb.py cg24860332 ...
        CPGS = [x.strip() for x in sys.argv[1:] if x.strip()]
    all_rows, meta = [], []
    for i, cpg in enumerate(CPGS, 1):
        cpg = cpg.strip()
        raw = os.path.join(RAWD, "%s.html" % cpg)
        rows, n = (None, 0)
        if os.path.exists(raw) and os.path.getsize(raw) > 5000:
            rows, n = parse_html(raw)
        attempt = 0
        while (rows is None or n == 0 or n <= 1) and attempt < 4:
            attempt += 1
            ok = curl_post(cpg, raw)
            rows, n = parse_html(raw)
            print("[%2d/%d] %s attempt=%d rc=%s rows=%s" % (i, len(CPGS), cpg, attempt, ok, n), flush=True)
            if rows is None:
                time.sleep(15 * attempt)
            elif n <= 1:
                break
        if rows is None:
            meta.append((cpg, "TRUNCATED", 0))
            print("  !! %s 响应被截断，重试仍失败" % cpg, flush=True)
            continue
        if n == 0:
            meta.append((cpg, "NO_RESULT", 0))
            continue
        # 写单 CpG
        outp = os.path.join(OUTD, "aries_mqtl_%s.tsv" % cpg)
        with io.open(outp, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\t".join(COLS) + "\n")
            for r in rows:
                if r[7] and r[7] != cpg:      # 服务端偶尔连带返回邻近 CpG，只留目标
                    continue
                fh.write("\t".join(r) + "\n")
        kept = sum(1 for r in rows if r[7] == cpg)
        meta.append((cpg, "OK", kept))
        for r in rows:
            if r[7] == cpg:
                all_rows.append(r)
        print("  -> %s 保留 %d 行" % (cpg, kept), flush=True)
        time.sleep(3)

    outp = os.path.join(OUTD, "aries_mqtl_all.tsv")
    if len(sys.argv) <= 1:                    # 分次抓取时不覆盖合并表
        with io.open(outp, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\t".join(COLS) + "\n")
            for r in all_rows:
                fh.write("\t".join(r) + "\n")
        print("\n写出 %s（%d 行）" % (outp, len(all_rows)), flush=True)
    sump = os.path.join(OUTD, "aries_fetch_status_%s.tsv" % ("batch" if len(sys.argv) > 1 else "all"))
    with io.open(sump, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("cpg\tstatus\tn_rows\n")
        for m in meta:
            fh.write("%s\t%s\t%d\n" % m)


if __name__ == "__main__":
    main()
