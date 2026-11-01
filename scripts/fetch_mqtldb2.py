#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fetch_mqtldb2.py —— mQTLdb(ARIES) 并发抓取器（带卡死检测）

为什么重写
----------
服务器端总带宽只有 ~20 KB/s，单连接串行抓 19 个 CpG 会被一个卡死的连接整块堵住
（实测 cg11348583 在 2.0 MB 处停了 5 分钟，后续 3 个已下载完的 CpG 排队不动）。
因此改为：4 线程 + `--speed-limit 2000 --speed-time 25` 自动掐断停滞连接 + 指数退避重试。

用法
----
python fetch_mqtldb2.py                # 抓默认 19 个 + 3 个 DEFA1 邻近代理 CpG
python fetch_mqtldb2.py cg24860332     # 只抓指定 CpG
"""
import os, re, sys, time, subprocess, io, threading, queue

BASE = "http://www.mqtldb.org/cgi-bin/search.cgi"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAWD = os.path.join(ROOT, "data", "mqtldb", "raw")
OUTD = os.path.join(ROOT, "results")
os.makedirs(RAWD, exist_ok=True)
os.makedirs(OUTD, exist_ok=True)

DEFAULT = [
    "cg07209298", "cg11348583", "cg13441112", "cg13682187",
    "cg00705004", "cg05969298", "cg07297068", "cg19304352",
    "cg24860332", "cg02787689", "cg05787038",          # DEFA1 簇内、紧邻 cg00705004 的代理 CpG
    "cg00539347", "cg27665648",
    "cg02846316", "cg03100639", "cg10451425", "cg04018474",
    "cg27535410", "cg23386659",
    "cg25599673", "cg17448192",
    "cg02671941",
]
COLS = ["timepoint", "snp", "snp_chr", "snp_pos", "a1", "a2", "maf",
        "cpg", "cpg_chr", "cpg_pos", "beta", "tstat", "effect_size", "pval", "trans"]

LOCK = threading.Lock()
LOG = []


def log(msg):
    with LOCK:
        LOG.append(msg)
        print(msg, flush=True)


def fetch(cpg, attempts=5):
    raw = os.path.join(RAWD, "%s.html" % cpg)
    for a in range(1, attempts + 1):
        tmp = raw + ".tmp"
        cmd = ["curl", "-s", "--ssl-no-revoke", "--max-time", "300",
               "--speed-limit", "2000", "--speed-time", "25",
               "--form-string", "batch=%s" % cpg,
               "-F", "database=1", "-F", "timepoint=0", "-F", "distance=1000000",
               BASE, "-o", tmp]
        subprocess.run(cmd, capture_output=True)
        ok = False
        if os.path.exists(tmp):
            with io.open(tmp, "r", encoding="utf-8", errors="replace") as fh:
                ok = "</table>" in fh.read()
            if ok:
                os.replace(tmp, raw)
                return True, a
            os.remove(tmp)
        log("  [retry] %s attempt=%d 未取完整，退避 %ds" % (cpg, a, 10 * a))
        time.sleep(10 * a)
    return False, attempts


def parse(path):
    with io.open(path, "r", encoding="utf-8", errors="replace") as fh:
        html = fh.read()
    rows = []
    for tr in re.findall(r"<tr>(.*?)</tr>", html, flags=re.S):
        tds = [re.sub(r"<[^>]+>", "", x).strip() for x in
               re.findall(r"<td[^>]*>(.*?)</td>", tr, flags=re.S)]
        if len(tds) >= 14:
            rows.append(tds[:15] if len(tds) >= 15 else tds + [""] * (15 - len(tds)))
    return rows


def worker(q, meta):
    while True:
        try:
            cpg = q.get_nowait()
        except queue.Empty:
            return
        raw = os.path.join(RAWD, "%s.html" % cpg)
        need = True
        if os.path.exists(raw) and os.path.getsize(raw) > 5000:
            with io.open(raw, "r", encoding="utf-8", errors="replace") as fh:
                need = "</table>" not in fh.read()
        if need:
            ok, n_att = fetch(cpg)
            if not ok:
                with LOCK:
                    meta.append((cpg, "FAILED", 0))
                log("  !! %s 抓取失败" % cpg)
                q.task_done()
                continue
        rows = [r for r in parse(raw) if r[7] == cpg]
        outp = os.path.join(OUTD, "aries_mqtl_%s.tsv" % cpg)
        with io.open(outp, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\t".join(COLS) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        with LOCK:
            meta.append((cpg, "OK", len(rows)))
        log("  -> %-12s %4d 行" % (cpg, len(rows)))
        q.task_done()


def main():
    cpgs = [x.strip() for x in sys.argv[1:]] or DEFAULT
    q = queue.Queue()
    for c in cpgs:
        q.put(c)
    meta = []
    ths = [threading.Thread(target=worker, args=(q, meta), daemon=True) for _ in range(4)]
    t0 = time.time()
    for t in ths:
        t.start()
    for t in ths:
        t.join()
    log("\n耗时 %.1f 分钟；结果：" % ((time.time() - t0) / 60))
    for c, s, n in sorted(meta):
        log("  %-12s %-8s %d" % (c, s, n))
    sump = os.path.join(OUTD, "aries_fetch_status_%s.tsv" % ("batch" if sys.argv[1:] else "all"))
    with io.open(sump, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("cpg\tstatus\tn_rows\n")
        for m in sorted(meta):
            fh.write("%s\t%s\t%d\n" % m)


if __name__ == "__main__":
    main()
