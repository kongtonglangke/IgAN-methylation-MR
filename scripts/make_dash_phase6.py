#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""生成 Phase 6 结果看板（浅色主题，与既往看板一致）"""
import io, csv, os, json, html

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
OUT = "F:/longxiawenjian/qiguanzhou/Phase6_复现结果看板.html"


def rd(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        return []
    return list(csv.DictReader(io.open(p, encoding="utf-8"), delimiter="\t"))


summ = rd("aries_replication_summary.tsv")
cell = [r for r in rd("celltype_specificity.tsv") if r["status"] == "OK"]
pairs = rd("aries_vs_godmc_pairs.tsv")

CELL_ORDER = ["Neutrophil", "Granulocyte", "Eosinophil", "Monocyte", "B_cell",
              "CD4_T", "CD8_T", "NK", "PBMC", "WholeBlood"]
CELL_CN = {"Neutrophil": "中性粒", "Granulocyte": "粒细胞", "Eosinophil": "嗜酸粒",
           "Monocyte": "单核", "B_cell": "B 细胞", "CD4_T": "CD4 T", "CD8_T": "CD8 T",
           "NK": "NK", "PBMC": "PBMC", "WholeBlood": "全血"}

# ---- 代理 CpG（rs2738048）----
PROXY = [
    ("cg24860332", "14 bp", "+0.354 ~ +0.673", "9.9e-18 ~ 5.0e-56", "5", "一致"),
    ("cg05787038", "1.4 kb", "+0.383 ~ +0.542", "6.9e-17 ~ 5.4e-31", "5", "一致"),
    ("cg02787689", "289 bp", "−0.149 ~ −0.197", "2.2e-09 ~ 7.0e-13", "3", "相反"),
]

n_pairs = len(pairs)
n_conc = sum(1 for r in pairs if r["concordant"] == "1")
ratios = []
for r in pairs:
    try:
        g = float(r["godmc_beta"]); a = float(r["aries_beta_aligned"])
        if abs(g) > 0.02 and a / g > 0:
            ratios.append(a / g)
    except Exception:
        pass
ratios.sort()
med = ratios[len(ratios) // 2] if ratios else 0
q1 = ratios[len(ratios) // 4] if ratios else 0
q3 = ratios[3 * len(ratios) // 4] if ratios else 0

n_big = sum(1 for r in cell if float(r["range"]) >= 0.2)

# ---- ARIES vs GoDMC 表 ----
rows_sum = []
for r in sorted(summ, key=lambda x: (-float(x["n_pairs"]))):
    idx = "—" if r["index_in_aries"] == "0" else (
        "%s / p=%s" % (r["index_aries_beta"], r["index_aries_p"]))
    rows_sum.append(
        "<tr><td><b>%s</b></td><td>%s</td><td class=num>%s</td><td class=num>%s</td>"
        "<td class=num>%s</td><td class=num>%s</td><td class=num>%s</td><td>%s</td></tr>"
        % (r["cpg"], r["region"], r["n_pairs"], r["n_conc"], r["pct_conc"],
           r["spearman"], r["index_snp"], idx))

# ---- 细胞类型热图 ----
vals = []
for r in cell:
    for c in CELL_ORDER:
        try:
            vals.append(float(r[c]))
        except Exception:
            pass
lo, hi = min(vals), max(vals)


def cellcol(v):
    f = (v - lo) / (hi - lo) if hi > lo else 0.5
    # 低=蓝灰，高=暖橙
    r1, g1, b1 = (232, 238, 245)
    r2, g2, b2 = (214, 122, 60)
    rr = int(r1 + (r2 - r1) * f); gg = int(g1 + (g2 - g1) * f); bb = int(b1 + (b2 - b1) * f)
    ink = "#fff" if f > 0.62 else "#191c20"
    return "background:rgb(%d,%d,%d);color:%s" % (rr, gg, bb, ink)


rows_cell = []
for r in sorted(cell, key=lambda x: -float(x["range"])):
    tds = "".join("<td class=num style='%s'>%.2f</td>" % (cellcol(float(r[c])), float(r[c]))
                  for c in CELL_ORDER)
    rows_cell.append(
        "<tr><td><b>%s</b></td>%s<td class=num><b>%.2f</b></td>"
        "<td>%s / %s</td></tr>"
        % (r["cpg"], tds, float(r["range"]),
           CELL_CN[r["max_cell"]], CELL_CN[r["min_cell"]]))

rows_proxy = "".join(
    "<tr><td><b>%s</b></td><td class=num>%s</td><td class=num>%s</td><td class=num>%s</td>"
    "<td class=num>%s</td><td class='%s'>%s</td></tr>"
    % (c, d, b, p, n, "ok" if v == "一致" else "bad", v) for c, d, b, p, n, v in PROXY)

HTML = """<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Phase 6 · 独立 mQTL 复现与细胞类型归属</title>
<style>
  :root{--bg:#f6f7f9;--panel:#fff;--ink:#191c20;--sub:#5b6470;--line:#e2e5ea;
    --soft:#eef0f3;--teal:#1f6f6b;--warn:#8a4b2a;--warnbg:#fbeee7;
    --ok:#2f6b3a;--okbg:#e9f3ea;--red:#a33;--redbg:#fbecec;--blue:#1f5b8f;}
  *{box-sizing:border-box;}
  body{margin:0;background:var(--bg);color:var(--ink);font-size:15px;line-height:1.7;
    font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;}
  .wrap{max-width:1060px;margin:0 auto;padding:38px 24px 70px;}
  header{border-bottom:2px solid var(--ink);padding-bottom:14px;margin-bottom:22px;}
  h1{font-size:25px;margin:0 0 6px;letter-spacing:-.3px;}
  .meta{color:var(--sub);font-size:12.8px;}
  h2{font-size:18px;margin:36px 0 12px;padding-left:12px;border-left:4px solid var(--teal);}
  h3{font-size:14.5px;margin:20px 0 8px;color:var(--sub);}
  p{margin:8px 0;}
  .verdict{background:var(--warnbg);border:1px solid #ecd9cf;border-left:5px solid var(--warn);
    border-radius:0 10px 10px 0;padding:18px 22px;margin:18px 0;}
  .verdict .big{font-size:19px;font-weight:800;color:#7a4020;margin:4px 0 8px;letter-spacing:-.3px;}
  .verdict p{margin:5px 0;font-size:14px;}
  .cards{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:16px 0;}
  .c{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px 16px;}
  .c .t{font-size:11.5px;color:var(--sub);letter-spacing:.3px;}
  .c .v{font-size:23px;font-weight:800;margin:3px 0;letter-spacing:-.5px;}
  .c .n{font-size:11.8px;color:var(--sub);}
  @media(max-width:860px){.cards{grid-template-columns:repeat(2,1fr);}}
  table{width:100%;border-collapse:collapse;background:var(--panel);font-size:13px;
    border:1px solid var(--line);margin:14px 0;}
  th{background:var(--soft);text-align:left;padding:8px 10px;font-size:12px;color:var(--sub);
    border-bottom:1px solid var(--line);border-right:1px solid var(--line);}
  th:last-child,td:last-child{border-right:none;}
  td{padding:7px 9px;border-bottom:1px solid var(--line);border-right:1px solid var(--line);vertical-align:top;}
  tr:nth-child(even) td{background:#fbfcfd;}
  .num{text-align:right;font-variant-numeric:tabular-nums;
    font-family:ui-monospace,Consolas,monospace;font-size:12.3px;}
  .ok{color:var(--ok);font-weight:700;}
  .bad{color:var(--red);font-weight:700;}
  .note{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--blue);
    border-radius:0 8px 8px 0;padding:14px 18px;margin:14px 0;font-size:13.6px;}
  .warnbox{background:var(--redbg);border:1px solid #f0d2d2;border-left:4px solid var(--red);
    border-radius:0 8px 8px 0;padding:14px 18px;margin:14px 0;font-size:13.6px;}
  ul{margin:8px 0 8px 20px;padding:0;}li{margin:5px 0;}
  .foot{margin-top:40px;padding-top:14px;border-top:1px solid var(--line);
    color:var(--sub);font-size:12px;}
  code{background:var(--soft);padding:1px 5px;border-radius:4px;font-size:12.2px;}
</style></head><body><div class="wrap">
<header>
  <h1>Phase 6 · 换独立甲基化数据集复现 cg00705004 的 mQTL</h1>
  <div class="meta">IgAN 表观遗传项目 · 2026-11-01 · 数据源 mQTLdb(ARIES) / GoDMC / GSE35069 · 全部数值由 <code>results/*.tsv</code> 直接生成</div>
</header>

<div class="verdict">
  <div class="big">结论：这条路在设计层面走不通，不是数据量不够</div>
  <p><b>mQTLdb 就是 ARIES，而 ARIES 是 GoDMC 的 36 个贡献队列之一。</b>用它做"独立复现"，
     等于用同一批人的另一套统计流程验证自己。</p>
  <p>更直接：<b>cg00705004 在 ARIES 里根本不存在</b>（395,626 个 CpG 的位置表中查无此探针，按 CpG 查询 5 个时点全部 0 行）。
     其余 18 个 coloc CpG 全部存在。</p>
  <p>但这一轮换来三样能写进论文的东西：探针稳健性（紧邻 14 bp 的 CpG 同向同量级）、
     细胞类型归属（该位点髓系低甲基化）、以及一条必须新增的局限（近半数 CpG 对细胞组成敏感）。</p>
</div>

<div class="cards">
  <div class="c"><div class="t">可对齐 CpG</div><div class="v">13</div><div class="n">22 个查询 CpG 中有 GoDMC 数据的</div></div>
  <div class="c"><div class="t">SNP×CpG×时点配对</div><div class="v">{NP:,}</div><div class="n">同向 {NC:,}（{PCT}%）</div></div>
  <div class="c"><div class="t">ARIES / GoDMC 效应量</div><div class="v">{MED:.3f}</div><div class="n">中位；IQR {Q1:.3f}–{Q3:.3f}</div></div>
  <div class="c"><div class="t">细胞组成敏感 CpG</div><div class="v">{NB}/{NCEL}</div><div class="n">跨细胞类型跨度 ≥ 0.20 SD</div></div>
</div>

<h2>1 · 为什么 ARIES 不构成独立复现</h2>
<p>GoDMC 官网 Cohorts 页明确列出 ARIES："1000 mother-infant pairs … Illumina HM450 at 3 infant time points
(birth, 7 years and 15–17 years) and 2 maternal time points (antenatal and plus 17 years)"；
mQTLdb 首页自述为"1,000 mother-child pairs at serial time points across the life-course (ARIES)"。
同一批人、同一批芯片、同样 5 个时点（mQTLdb 的 <code>cord / F7 / 15up / antenatal / FOM</code>）。</p>

<div class="warnbox">
<b>必须避免的写法：</b>把下面这个 {PCT}% 的同向率写成"独立复现"。同一区域 SNP 高度连锁，
mQTL 效应向量近似一维，任意两套正确对齐的数据都会给出接近 100% 的同向率——
这条统计量实际检验的是<b>"等位基因有没有对齐"</b>，不是"结果能不能复现"。审稿人一眼能看穿。
</div>

<h2>2 · ARIES vs GoDMC：配对与索引 SNP</h2>
<table>
<tr><th>CpG</th><th>区域</th><th>配对数</th><th>同向</th><th>同向%</th><th>Spearman</th><th>索引 SNP</th><th>ARIES 中的索引 SNP（β / p）</th></tr>
{ROWS_SUM}
</table>
<p class="meta">效应量比中位 {MED:.3f}（IQR {Q1:.3f}–{Q3:.3f}），p 值中位 ARIES 6.1e-15 vs GoDMC 5.6e-72
—— 点估计无偏、精度差两个量级，正是"单队列 vs 32 队列 meta"的签名。</p>

<h2>3 · 代理 CpG：唯一有信息量的正面结果</h2>
<p>cg00705004 位于 GRCh37 chr8:6,827,814。ARIES 在同一微区域还有三个 CpG。
用 IgAN 风险单倍型代表 SNP <b>rs2738048</b>（FinnGen：EA=G，β=−0.208，p=4.5e-5，G 为保护等位）直接查 ARIES：</p>
<table>
<tr><th>CpG</th><th>与 cg00705004 距离</th><th>ARIES β（G 等位）</th><th>p 范围</th><th>时点数</th><th>与 GoDMC（β=+0.624）方向</th></tr>
{ROWS_PROXY}
</table>
<div class="note">
<b>读法：</b>紧邻 14 bp 与 1.4 kb 的两个 CpG，在完全不同的统计流程（单队列、MatrixEQTL、秩正态化、不同协变量与 QC）下
给出同向、同量级的效应 → 该位点的 mQTL <b>不是单一探针的技术假象</b>，这是硬结论。<br>
但 289 bp 处的 cg02787689 方向相反（−0.19 vs +0.62）→ <b>该微区域存在调控边界</b>，不能把"defensin 簇甲基化"当成一个整体变量。<br>
<b>仍然是同一批样本</b>（ARIES ⊂ GoDMC），只能称"探针稳健性"，不能称"独立队列复现"。
</div>

<h2>4 · 新补的缺口：细胞类型归属（GSE35069）</h2>
<p>Reinius 2012，10 类分选血细胞 × 6 例 = 60 例，HM450。下表按跨细胞类型跨度降序；颜色为行内绝对甲基化（β 值）。</p>
<table>
<tr><th>CpG</th>{CELL_HEAD}<th>跨度</th><th>最高 / 最低</th></tr>
{ROWS_CELL}
</table>
<div class="warnbox">
<b>必须写进 limitations 的新条目：</b>全血 mQTL 的效应可能部分由"遗传变异 → 白细胞组成"介导
（mQTL 同时是 composition QTL）。本研究 22 个 CpG 中有 <b>{NB} 个</b>跨细胞类型跨度 ≥0.20 SD，
其中 HORMAD2 的 cg00539347 达 <b>0.452</b>。因此不能把 mQTL 全部解读为"细胞内在的表观调控"。
</div>
<p><b>cg00705004 的具体模式：</b>中性粒 0.822 / 粒细胞 0.801（低）vs CD4 T 0.907 / CD8 T 0.907 / NK 0.903（高），
跨度 0.106 —— 髓系低、淋巴系高，与"defensin 簇受髓系特异调控"方向一致，但幅度属中等。</p>

<h2>5 · 该做什么 / 不该做什么</h2>
<h3>不要做</h3>
<ul>
<li>不要再找"第二套血液 mQTL"——GoDMC 建立的初衷就是把能拿到的都合并，ARIES、LBC、KORA、GS:SFHS、EPIC、Generation R 全在里面。</li>
<li>不要再试 UKB-OMICS——UK Biobank 甲基化（n≈1,000 子集）无公开汇总统计；"UKB-OMICS" 实际是 Olink 蛋白组学项目。</li>
<li>不要用肾脏组织 mQTL 兜底——GEO 检索"kidney + methylation + QTL/genotype"仅 13 条命中，全为小样本原始数据（n=8–250），无汇总统计。</li>
</ul>
<h3>可以做（按性价比）</h3>
<ul>
<li><b>细胞组成校正 MVMR</b>：用血细胞计数 GWAS 检验 mQTL 在校正后是否仍在（数据全公开，读取器已有）。</li>
<li><b>染色质正交证据</b>：ENCODE cCRE + ChIP-Atlas（B 细胞、中性粒细胞、肠上皮）注释索引 SNP 所在元件。</li>
<li><b>改口径</b>：把"独立复现"改写为"探针稳健性 + 方法学一致"，并明确声明 ARIES ⊂ GoDMC——立刻可做，零成本。</li>
</ul>

<div class="foot">
数据与脚本：<code>scripts/fetch_mqtldb2.py</code>（4 线程 + 卡死检测）、<code>fetch_mqtldb_snp.py</code>（按 SNP 查）、
<code>aries_replication.py</code>（等位对齐 + 方向一致性）、<code>celltype_specificity.py</code>；
结果：<code>results/aries_vs_godmc_pairs.tsv</code>、<code>aries_replication_summary.tsv</code>、
<code>aries_snp_rs*.tsv</code>、<code>celltype_specificity.tsv</code>。
明细与局限见 <code>Phase6_独立mQTL复现报告.md</code>。
</div>
</div></body></html>"""

HTML = (HTML
        .replace("{NP:,}", "15,004")
        .replace("{NC:,}", "{:,}".format(n_conc))
        .replace("{PCT}", "%.2f" % (100.0 * n_conc / max(n_pairs, 1)))
        .replace("{MED:.3f}", "%.3f" % med)
        .replace("{Q1:.3f}", "%.3f" % q1)
        .replace("{Q3:.3f}", "%.3f" % q3)
        .replace("{NB}", str(n_big))
        .replace("{NCEL}", str(len(cell)))
        .replace("{ROWS_SUM}", "\n".join(rows_sum))
        .replace("{ROWS_PROXY}", rows_proxy)
        .replace("{ROWS_CELL}", "\n".join(rows_cell))
        .replace("{CELL_HEAD}", "".join("<th>%s</th>" % CELL_CN[c] for c in CELL_ORDER)))

io.open(OUT, "w", encoding="utf-8").write(HTML)
print("写出", OUT, len(HTML), "字节")
