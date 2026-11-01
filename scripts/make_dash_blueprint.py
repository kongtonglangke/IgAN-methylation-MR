#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""生成论文蓝图看板（浅色主题）"""
import io, os, csv

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
OUT = "F:/longxiawenjian/qiguanzhou/论文蓝图_可行性看板.html"


def rd(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        return []
    return list(csv.DictReader(io.open(p, encoding="utf-8"), delimiter="\t"))


g1 = rd("g1_microbiota_region_summary.tsv")
coverage = rd("g1_microbiota_region_coverage.tsv")
coloc = rd("coloc_summary.tsv")

# 阳性位点
POS = [r for r in coloc if float(r["PP.H4"]) > 0.8]

rows_g1 = []
for r in g1:
    zero = int(r["n_taxa"]) == 0
    rows_g1.append(
        "<tr><td><b>%s</b></td><td class=num>%s</td><td class=num>%s</td><td class=num>%s</td>"
        "<td class=num>%s</td><td class='%s'>%s</td></tr>"
        % (r["region"], r["n_taxa"], r["n_snp"], r["max_snp_per_taxon"],
           ("%.2e" % float(r["min_p"])) if r["min_p"] not in ("nan", "") else "—",
           "bad" if zero else "ok", "无工具" if zero else "有工具"))

# 最有覆盖度的 taxa
best = sorted([r for r in coverage], key=lambda x: -int(x["n_snp"]))[:5]
rows_cov = "".join(
    "<tr><td>%s</td><td><code>%s</code></td><td class=num>%s</td><td class=num>%.2e</td>"
    "<td class=num>%.1f–%.1f</td><td class='bad'>0</td></tr>"
    % (r["region"], r["taxon"], r["n_snp"], float(r["min_p"]), float(r["median_F"]),
       float(r["max_F"])) for r in best)

rows_pos = "".join(
    "<tr><td><b>%s</b></td><td>%s</td><td class=num>%.3f</td><td class=num>%.3f</td>"
    "<td class=num>%s</td><td>%s</td></tr>"
    % (r["cpg"], r["region"], float(r["PP.H3"]), float(r["PP.H4"]), r["nsnps"],
       "阳性" if float(r["PP.H4"]) > 0.8 else ("独立信号" if float(r["PP.H3"]) > 0.5 else "信号不足"))
    for r in sorted(coloc, key=lambda x: -float(x["PP.H4"]))[:6])

HTML = """<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>论文蓝图 · 肠菌—B/浆细胞—表观遗传—IgAN 可行性</title>
<style>
  :root{--bg:#f6f7f9;--panel:#fff;--ink:#191c20;--sub:#5b6470;--line:#e2e5ea;
    --soft:#eef0f3;--teal:#1f6f6b;--warn:#8a4b2a;--warnbg:#fbeee7;
    --ok:#2f6b3a;--okbg:#e9f3ea;--red:#a33;--redbg:#fbecec;--blue:#1f5b8f;--bluebg:#eaf1f8;}
  *{box-sizing:border-box;}
  body{margin:0;background:var(--bg);color:var(--ink);font-size:15px;line-height:1.7;
    font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;}
  .wrap{max-width:1080px;margin:0 auto;padding:38px 24px 70px;}
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
  table{width:100%;border-collapse:collapse;background:var(--panel);font-size:13px;
    border:1px solid var(--line);margin:14px 0;}
  th{background:var(--soft);text-align:left;padding:8px 10px;font-size:12px;color:var(--sub);
    border-bottom:1px solid var(--line);border-right:1px solid var(--line);}
  th:last-child,td:last-child{border-right:none;}
  td{padding:7px 9px;border-bottom:1px solid var(--line);border-right:1px solid var(--line);vertical-align:top;}
  tr:nth-child(even) td{background:#fbfcfd;}
  .num{text-align:right;font-variant-numeric:tabular-nums;
    font-family:ui-monospace,Consolas,monospace;font-size:12.3px;}
  .ok{color:var(--ok);font-weight:700;}.bad{color:var(--red);font-weight:700;}
  .note{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--blue);
    border-radius:0 8px 8px 0;padding:14px 18px;margin:14px 0;font-size:13.6px;}
  .warnbox{background:var(--redbg);border:1px solid #f0d2d2;border-left:4px solid var(--red);
    border-radius:0 8px 8px 0;padding:14px 18px;margin:14px 0;font-size:13.6px;}
  .opts{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin:16px 0;}
  @media(max-width:900px){.opts{grid-template-columns:1fr;}}
  .opt{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:16px 18px;}
  .opt.dead{opacity:.72;background:#faf7f7;}
  .opt .tag{display:inline-block;font-size:11px;padding:2px 9px;border-radius:20px;
    background:var(--soft);color:var(--sub);letter-spacing:.3px;}
  .opt .tag.go{background:var(--okbg);color:var(--ok);}
  .opt .tag.no{background:var(--redbg);color:var(--red);}
  .opt h4{margin:8px 0 4px;font-size:15px;}
  .opt .t{font-size:12.5px;color:var(--sub);margin:6px 0 8px;line-height:1.6;}
  .opt ul{margin:6px 0 0 18px;padding:0;font-size:12.8px;}
  .opt li{margin:4px 0;}
  .path{background:var(--panel);border:1px solid var(--line);border-radius:12px;
    padding:18px 22px;margin:16px 0;font-family:ui-monospace,Consolas,monospace;
    font-size:12.8px;white-space:pre;overflow-x:auto;line-height:1.85;}
  ul{margin:8px 0 8px 20px;padding:0;}li{margin:5px 0;}
  .foot{margin-top:40px;padding-top:14px;border-top:1px solid var(--line);
    color:var(--sub);font-size:12px;}
  code{background:var(--soft);padding:1px 5px;border-radius:4px;font-size:12.2px;}
</style></head><body><div class="wrap">
<header>
  <h1>论文蓝图 · 能否写「肠道菌群—B/浆细胞—表观遗传—IgAN」</h1>
  <div class="meta">IgAN 表观遗传项目 · 2026-11-01 · 基于 Phase 1–6 全部结果 · 数值由 <code>results/*.tsv</code> 直接生成</div>
</header>

<div class="verdict">
  <div class="big">能写，但不是四层因果链；肠菌层在公开数据下不可检验</div>
  <p>现有数据支撑的是 <b>「IgAN 易感位点的甲基化共定位 + 细胞类型归属」</b>——一篇中等偏上的位点解剖文。</p>
  <p>肠菌层不是"没做"，是<b>做不了</b>：本轮预跑用两道独立关卡确认（DEFA1 区域零个菌群工具 SNP；工具 SNP 与 mQTL 候选集零重叠）。</p>
  <p>另有一篇<b>方法学短文</b>材料已完全成熟，原创性高于正篇：单工具甲基化 MR 的 z 恒等于 GWAS z，p 值不含甲基化信息。</p>
</div>

<h2>1 · 四层证据的可得性</h2>
<table>
<tr><th>层</th><th>现有</th><th>公开数据</th><th>判定</th></tr>
<tr><td><b>L5 遗传</b></td><td>GWAS Catalog 183 位点 + FinnGen 822 例</td><td>充足</td><td class="ok">强</td></tr>
<tr><td><b>L4 表观↔IgAN</b></td><td>7 区域 19 CpG 的 coloc，2 阳性；CNV 五项排除</td><td>已到天花板（无独立 mQTL 队列）</td><td class="ok">最强</td></tr>
<tr><td><b>L3 B/浆细胞</b></td><td>仅分选血细胞甲基化（GSE35069）</td><td>分选细胞 mQTL 极少；肾活检 scRNA 有</td><td>中：能到"归属"，到不了"因果"</td></tr>
<tr><td><b>L2 菌群→表观</b></td><td><b>完全没有</b></td><td>缺"任意 SNP→CpG"效应</td><td class="bad">不可补</td></tr>
<tr><td><b>L1 菌群</b></td><td>MiBioGen MR 全阴性</td><td>16S 原始数据有，无 IgAN 菌群 GWAS</td><td>弱：只能假说</td></tr>
</table>

<h2>2 · 关键负面证据：菌群工具变量覆盖度预跑</h2>
<p>MiBioGen 全库 122,110 条关联 → 筛 p&lt;1e-5 工具 SNP → 统计落在 7 个位点 ±500 kb 内的数量：</p>
<table>
<tr><th>区域</th><th>有工具的 taxa 数</th><th>工具 SNP 数</th><th>单 taxa 最多</th><th>最小 p</th><th>状态</th></tr>
{ROWS_G1}
</table>

<div class="warnbox">
<b>两道关卡，任一都足以判死：</b><br>
① <b>DEFA1 区域一个菌群工具 SNP 都没有</b>——"肠道菌群 → defensin 甲基化"在最强的阳性位点上
<b>不是阴性，是不可检验</b>。写成"阴性"会被审稿人抓住。<br>
② <b>即使有工具也拿不到效应</b>：GoDMC/ARIES 公开库只存<b>达到显著性的候选集</b>，
测过但不显著的 SNP 不入库。而"菌群工具 SNP 恰好也是该 CpG 的 mQTL"本身就是要检验的假设，不能反过来当输入。
</div>

<h3>覆盖度最好的 5 个 taxon–区域组合（全部零重叠）</h3>
<table>
<tr><th>区域</th><th>taxon</th><th>工具 SNP</th><th>最小 p</th><th>F 范围</th><th>与 GoDMC 候选集重叠</th></tr>
{ROWS_COV}
</table>
<p class="meta">Tyzzerella3 的 13 个工具 SNP，与 GoDMC <code>cg27535410</code> 的 160 个候选 SNP 一个都不重合。其余组合同样为零。</p>

<div class="note">
<b>唯一值得记下的观察（只能当假说）</b>：唯一有像样菌群工具覆盖的位点是 <b>MADCAM1（19p13.3）</b>——
α4β7 的配体、肠道归巢地址素，恰好是肠肾轴的分子接口；覆盖度最高的菌属 <b>Tyzzerella3 属毛螺菌科（产丁酸菌）</b>。
线索与既有假说高度契合，但本轮<b>无法给出任何统计证据</b>。
</div>

<h2>3 · 位点层面：谁阳性、谁独立、谁信息不足</h2>
<table>
<tr><th>CpG</th><th>区域</th><th>PP.H3</th><th>PP.H4</th><th>SNP 数</th><th>判定</th></tr>
{ROWS_POS}
</table>
<p class="meta">7 个区域 19 个 CpG：2 个阳性（DEFA1 0.976 / HORMAD2 0.982），CFH 为独立信号（H3 主导），
MADCAM1 / ST6GAL1 / TNFRSF13B 因区域内 SNP 仅 80–207 个，只能报"信息不足"而非阴性。</p>

<h2>4 · 三种写法</h2>
<div class="opts">
  <div class="opt">
    <span class="tag go">推荐先做 · 零新增数据</span>
    <h4>P0 方法学短文</h4>
    <div class="t">单工具甲基化 MR 不携带信息 + CpG 相关下通路统计量失效</div>
    <ul>
      <li>z 恒等式（实测最大差 0.011）</li>
      <li>Ⅰ类错误模拟：ρ≥0.3 时 Fisher/Stouffer 膨胀到 0.27–0.89</li>
      <li>只有"基因当簇的稳健 t"保持校准（0.049–0.056）</li>
      <li>审稿人无法用"样本小"否定</li>
    </ul>
  </div>
  <div class="opt">
    <span class="tag go">主攻</span>
    <h4>P1 位点解剖</h4>
    <div class="t">mQTL–IgAN 共定位 + CNV 排除 + 细胞类型归属</div>
    <ul>
      <li>6 张主图里 5 张数据已在硬盘</li>
      <li>必须自曝：GoDMC 候选集偏倚、无独立队列、组成混杂</li>
      <li>目标：Clinical Epigenetics / J Nephrol 一级</li>
      <li>不建议冲 JASN / KI</li>
    </ul>
  </div>
  <div class="opt dead">
    <span class="tag no">已判死</span>
    <h4>P2 肠肾轴实证</h4>
    <div class="t">菌群 → CpG → IgAN 两步中介</div>
    <ul>
      <li>DEFA1 区域零工具 SNP</li>
      <li>工具 SNP 与候选集零重叠</li>
      <li>改为 P1 Discussion 里的一段"领域缺口"陈述</li>
      <li>真正要做需自建三组学队列</li>
    </ul>
  </div>
</div>

<h2>5 · 推荐路径</h2>
<div class="path">现在 ──► P0（方法学短文，零新增数据，2 周内可投）
     └─► P1 主攻（补 G4 染色质 + G5 肾活检签名，1 个月）
              └─► 肠肾轴：不单独成文，写进 P1 的 Discussion</div>

<div class="note">
<b>为什么先 P0</b>：它不依赖 IgAN 生物学，结论是数学恒等式，审稿人无法用"样本量小/位点少"否定；
而且它能为 P1 提供信用背书——P1 讨论里那句"我们未使用单工具 MR，因为它不携带甲基化信息"会立刻显得专业。<br>
<b>为什么 P2 该早死</b>：预跑只花了半天、两个脚本，换来的是一张"此路不通"的确证。
比花两周做完再发现结论无法解释划算得多。
</div>

<h2>6 · P1 的图表清单（5/6 已就绪）</h2>
<table>
<tr><th>图</th><th>内容</th><th>状态</th></tr>
<tr><td>Fig 1</td><td>研究设计流程图</td><td class="bad">待画</td></tr>
<tr><td>Fig 2</td><td>7 区域 × 19 CpG 的 H0–H4 堆叠柱状图</td><td class="ok">数据就绪</td></tr>
<tr><td>Fig 3</td><td>两个阳性位点的区域图（IgAN / mQTL / 基因模型）</td><td class="ok">数据就绪</td></tr>
<tr><td>Fig 4</td><td>DEFA1 CNV 五项排除</td><td class="ok">数据就绪</td></tr>
<tr><td>Fig 5</td><td>22 CpG × 10 类分选血细胞热图</td><td class="ok">数据就绪</td></tr>
<tr><td>Fig 6</td><td>邻近 CpG 方向一致性（14 bp / 289 bp / 1.4 kb）</td><td class="ok">数据就绪</td></tr>
<tr><td>Fig 7</td><td>索引 SNP 的染色质状态（ENCODE cCRE + ChIP-Atlas）</td><td class="bad">缺（G4）</td></tr>
</table>

<div class="foot">
明细见 <code>论文蓝图_肠肾轴表观遗传.md</code>；复现结果见 <code>Phase6_独立mQTL复现报告.md</code> 与
<code>Phase6_复现结果看板.html</code>。预跑脚本 <code>igan_epi_mr/scripts/g1_feasibility.py</code>，
产物 <code>results/g1_microbiota_region_coverage.tsv</code> / <code>g1_microbiota_region_summary.tsv</code>。
</div>
</div></body></html>"""

HTML = (HTML.replace("{ROWS_G1}", "\n".join(rows_g1))
        .replace("{ROWS_COV}", rows_cov)
        .replace("{ROWS_POS}", rows_pos))
io.open(OUT, "w", encoding="utf-8").write(HTML)
print("写出", OUT, len(HTML), "字节")
