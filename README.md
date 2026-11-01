# Code and derived data

**When single-instrument methylation MR is uninformative: a methodological autopsy and colocalization reanalysis at IgA nephropathy susceptibility loci**

Yan Wang — Department of Nephrology, Affiliated Hospital of Gansu Medical College, Pingliang, Gansu 744000, China

This repository contains every analysis script, every derived result table behind the figures and
tables of the paper, and the figure source files. Raw third-party datasets are **not** redistributed
here; they are publicly available from their original sources, listed in `DATA_SOURCES.md`.

## Layout

| Path | Contents |
|---|---|
| `scripts/` | 57 analysis and figure scripts (Python 3 and R) |
| `results/` | 60 derived tables produced by those scripts |
| `figures/` | Fig 1–8 as PDF (vector) and PNG (400 dpi) |
| `tables/` | Supplementary Table S1 |

## Figure map

| Figure | Script |
|---|---|
| Fig 1 | `scripts/fig1_design.R` |
| Fig 2, Fig 4 | `scripts/fig2_fig3.R` |
| Fig 3 | `scripts/fig3_forest.R` |
| Fig 5 | `scripts/fig5_defa1_cnv.R` |
| Fig 6 | `scripts/fig6_celltype.R` |
| Fig 7 | `scripts/fig7_chromatin.R` |
| Fig 8 | `scripts/fig8_gdia1.R` |

Fig 2 and Fig 4 share one script because they derive from the same colocalization assembly step.

## Table map

| Table | Source |
|---|---|
| Table 1 | `scripts/sim_type1.py` (simulation) |
| Table 2 | `scripts/coloc_igan_mqtl.R` with `scripts/coloc_robustness.R` |
| Table 3 | `scripts/fig5_defa1_cnv.R` (segmental duplication, structural variation, colocalization robustness) and `scripts/cnv_signature.py` (dose signature) |
| Table 4 | Literature summary, no computation |
| Table S1 | `scripts/g1_feasibility.py` |

## Environment

- R 4.3.3 (2024-02-29), with ggplot2 4.0.3, coloc 5.2.3
- Python 3.13
- Remote tabix reads used a custom HTTP Range reader (`scripts/tabix_http.py`) because the mQTL
  resources are served as remote tabix files rather than bulk downloads.

## Reproducing the simulation in Table 1

```
python scripts/sim_type1.py
```

The script fixes the random seed (20261030) and runs all eight rho x density cells in sequence.
Running a single cell in isolation will not reproduce the published percentages exactly, because
the random-number stream then differs; see Supplementary methods S1.5.

## Script index

Descriptions are quoted verbatim from each script's header, so some remain in the author's working
language.

| Script | Header description (quoted verbatim) |
|---|---|
| `analyze_stage1b.R` | Stage 1b summary: CpG methylation -> IgA nephropathy (Sakaue 2021, GCST90018866) |
| `analyze_stage1b.py` | Summarise stage-1b (CpG methylation -> IgA nephropathy). |
| `analyze_stage2_pct.R` | Stage 2b: the calibrated pathway-level statistic applied to the real data. |
| `aries_replication.py` | aries_replication.py —— ARIES(mQTLdb) vs GoDMC 的独立队列方向一致性检验 |
| `celltype_specificity.py` | celltype_specificity.py —— 用 GSE35069（Reinius 2012, 纯化血细胞 450K, n=60） |
| `cnv_confound.py` | DEFA1 CNV-confounding diagnostics. |
| `cnv_signature.py` | Is the DEFA1 mQTL a copy-number signature or a local regulatory effect? |
| `coloc_cnv_sensitivity.R` | DEFA1 colocalisation: does the result survive removal of the duplicated |
| `coloc_igan_mqtl.R` | Colocalisation of IgAN risk loci with target-gene methylation |
| `coloc_prior_sensitivity.R` | Prior sensitivity of PP.H4 for the two colocalized CpGs |
| `coloc_robustness.R` | Robustness of the two positive colocalisations. |
| `compute_f.py` | Compute per-CpG instrument strength from GoDMC cis-mQTL. |
| `diag_igan_density.py` | Diagnose whether GCST90018866 (Sakaue 2021 IgAN) harmonised file is |
| `dl.py` | Multi-threaded range downloader for EBI GWAS summary statistics. |
| `dl_single.py` | Single-stream resumable downloader. Safe to re-run: always resumes from the |
| `fetch_chunked.py` | Fetch GoDMC mQTL in small chunks to avoid the server's large-query 403 limit. |
| `fetch_coloc_finngen.py` | Fetch FinnGen R12 IgA nephropathy summary statistics for the coloc regions. |
| `fetch_coloc_igan.py` | Fetch IgAN summary statistics around each colocalisation region. |
| `fetch_coloc_mqtl.py` | Fetch full cis-mQTL records for the CpGs that colocalise with IgAN risk loci. |
| `fetch_igan_loci.py` | Fetch every GWAS Catalog association for IgA nephropathy (MONDO_0005342) |
| `fetch_igan_lookup.py` | Fetch IgAN (Sakaue 2021, GCST90018866) summary statistics for the cis-mQTL |
| `fetch_mqtl.py` | X/Y 跳过（GoDMC 主要常染色体） |
| `fetch_mqtl2.py` | Fetch GoDMC cis-mQTL for candidate regions (GRCh37). |
| `fetch_mqtldb.py` | fetch_mqtldb.py —— 从 mQTLdb (ARIES, Gaunt 2016 Genome Biol) 取指定 CpG 的 mQTL 结果 |
| `fetch_mqtldb2.py` | fetch_mqtldb2.py —— mQTLdb(ARIES) 并发抓取器（带卡死检测） |
| `fetch_mqtldb_snp.py` | fetch_mqtldb_snp.py —— 按 SNP 查 mQTLdb(ARIES)，返回该 SNP 的全部 cis-mQTL CpG |
| `fetch_rest.py` | Sequential (throttled) fetch of remaining GoDMC regions. |
| `fig1_design.R` | fig1_design.R —— Fig 1：研究设计与数据流（示意图，无数据） |
| `fig2_fig3.R` | fig2_fig3.R —— 生成 Fig 2（恒等式散点 + Ⅰ类错误热图）与 Fig 3（coloc 后验堆叠图） |
| `fig3_forest.R` | fig3_forest.R —— Fig 3：16 个 IgAN 风险等位 → 血液甲基化（随机效应）森林图 |
| `fig5_defa1_cnv.R` | fig5_defa1_cnv.R —— Fig 5：DEFA1 位点的拷贝数 / 重复 / 探针特异性检查（hg19） |
| `fig6_celltype.R` | fig6_celltype.R —— Fig 6：22 个 CpG 在分选血细胞中的甲基化（GSE35069） |
| `fig7_chromatin.R` | fig7_chromatin.R —— Fig 7：两个共定位 CpG 的染色质注释（hg38） |
| `fig8_gdia1.R` | fig8_gdia1.R —— Fig 8：O-糖基化通路（C1GALT1 / C1GALT1C1）的甲基化 QTL |
| `finish_dl.py` | Append the missing tail of specific parts, then merge and verify gzip integrity. |
| `fix_parts.py` | Re-download only the parts whose size does not match, then merge and verify gzip. |
| `g1_feasibility.py` | g1_feasibility.py —— 「菌群 → CpG」效应矩阵的可行性预跑 |
| `g4_atac_direct.py` | G4: 直读 ChIP-Atlas ATAC 的指定区域 |
| `g4_atac_summary.py` | G4: ATAC peak 的细胞类群归属与最近距离 |
| `g4_ccre.py` | G4: ENCODE cCRE 注释 + 区域调控密度对比 |
| `g4_ccre_baseline.py` | G4 补充: 随机窗口 cCRE 密度基线 |
| `g4_chipatlas_atac.py` | G4: ChIP-Atlas hg38 ATAC-seq (ATC.ALL.05.AllAg.AllCell) 区域随机读取 |
| `get_coords.py` | — |
| `install_anno.R` | — |
| `liftover_chain.py` | Offline GRCh37 -> GRCh38 point liftOver from the UCSC chain file. |
| `make_dash_blueprint.py` | 生成论文蓝图看板（浅色主题） |
| `make_dash_phase6.py` | 生成 Phase 6 结果看板（浅色主题，与既往看板一致） |
| `match_igan_mqtl.py` | Option A core test: are the known IgA nephropathy risk variants themselves |
| `mr_cpg_bcell.py` | Stage 1 gate test: cis-mQTL of candidate CpGs -> outcome (B-cell trait or IgAN). |
| `mr_cpg_igan.py` | Stage 1b: CpG methylation -> IgA nephropathy (IgA 肾病) causal test. |
| `optionA_direction.py` | Option A, step 3: align the IgAN risk allele with the GoDMC effect allele and |
| `resume_merge.py` | Resume the final partial chunk and merge all parts into the output file. |
| `retry_coloc_mqtl.py` | Retry the CpGs that the GoDMC API refused with 403 during the bulk pull. |
| `scan_orru.py` | — |
| `scan_orru2.py` | Scan Orru 2020 immune-cell trait metadata on EBI FTP. |
| `sim_type1.py` | Stage 2a: type-I error simulation for candidate pathway-level statistics. |
| `tabix_http.py` | Minimal tabix/BGZF random-access reader over plain HTTPS Range requests. |

## Citation

If you use this code, please cite the paper (see `CITATION.cff`) and the Zenodo archive.
