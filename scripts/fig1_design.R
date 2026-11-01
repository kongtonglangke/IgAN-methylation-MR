#!/usr/bin/env Rscript
# fig1_design.R —— Fig 1：研究设计与数据流（示意图，无数据）
#
# 三栏：数据源 → 分析步骤（编号 1-8，对应正文 §3.1-3.8）→ 关键发现
# 输出：figs/Fig1.pdf、figs/Fig1.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

ROOT <- Sys.getenv("IGAN_MR_ROOT", unset = ".")
setwd(ROOT)
dir.create("figs", showWarnings = FALSE)

steps <- data.frame(
  i = 1:8,
  input = c(
    "GoDMC + ARIES mQTL\nFinnGen R12 IgAN GWAS",
    "Within-gene CpG\ncorrelation rho",
    "GWAS Catalog\n(183 unique IgAN loci)",
    "FinnGen + mQTL\n(shared SNP sets)",
    "UCSC segmental duplications\n1000G SV; 50-nt probe k-mer",
    "GSE35069; ENCODE cCRE\nChIP-Atlas ATAC",
    "Gd-IgA1 GWAS\n(Kiryluk et al. 2017)",
    "MiBioGen\n(456 instruments, 211 taxa)"
  ),
  analysis = c(
    "1. Single-CpG MR (1,699 pairs)\nWald ratio vs outcome z",
    "2. Type I error simulation\n6 strategies, 2,000 replicates/cell",
    "3. Instrument-first pairing\n32 candidate pairs",
    "4. Colocalization (coloc.abf)\n19 CpGs, 7 regions",
    "5. CNV / repeat / probe checks\n(DEFA1 locus)",
    "6. Cell-type + chromatin\nannotation",
    "7. C1GALT1 pathway mQTL\n(5 developmental windows)",
    "8. Microbiota instrument\ncoverage assessment"
  ),
  finding = c(
    "Identity: max |dz| = 7.4e-4\nNo methylation information",
    "Naive aggregation: 0.885\nGene-cluster t: 0.056",
    "16 pairs retained\nCFH fixed-effect artefact out",
    "2 / 19 colocalized\ncg00705004, cg00539347",
    "DEFA1 signal\nis not an artefact",
    "dELS + 1,752 blood ATAC\nvs no cCRE + 0 ATAC",
    "Robust mQTL, but\nmethylation != mediation",
    "0 / 456 pairable\nUntestable, not negative"
  ),
  bad = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE)
)
steps$y <- 8.55 - (steps$i - 1) * 1.06

XL <- 0.05; XC <- 3.15; XR <- 6.85   # 三栏中心

p <- ggplot(steps) +
  # 栏标题
  annotate("text", x = 1.55, y = 9.35, label = "Data sources", size = 3.2,
           fontface = "bold", colour = "grey20") +
  annotate("text", x = 4.95, y = 9.35, label = "Analysis (Results section)", size = 3.2,
           fontface = "bold", colour = "grey20") +
  annotate("text", x = 8.45, y = 9.35, label = "Finding", size = 3.2,
           fontface = "bold", colour = "grey20") +
  # 输入 → 分析 箭头
  geom_segment(aes(x = 2.72, xend = 3.02, y = y, yend = y),
               arrow = arrow(length = unit(1.6, "mm")), colour = "grey55",
               linewidth = 0.45) +
  # 分析 → 发现 箭头
  geom_segment(aes(x = 6.42, xend = 6.72, y = y, yend = y),
               arrow = arrow(length = unit(1.6, "mm")), colour = "grey55",
               linewidth = 0.45) +
  # 分析框
  geom_rect(aes(xmin = 3.05, xmax = 6.38, ymin = y - 0.40, ymax = y + 0.40),
            fill = "#EDF2F8", colour = "#4A7BB7", linewidth = 0.45) +
  # 步骤间竖向箭头
  geom_segment(data = steps[-8, ],
               aes(x = 4.72, xend = 4.72, y = y - 0.42, yend = y - 0.64),
               arrow = arrow(length = unit(1.6, "mm")), colour = "#4A7BB7",
               linewidth = 0.45) +
  # 文字
  geom_text(aes(x = 1.55, y = y, label = input), size = 2.35, colour = "grey35",
            lineheight = 0.95) +
  geom_text(aes(x = 4.72, y = y, label = analysis), size = 2.5, colour = "grey10",
            lineheight = 0.98, fontface = "bold") +
  geom_text(aes(x = 8.45, y = y, label = finding,
                colour = ifelse(bad, "warn", "ok")), size = 2.4, lineheight = 0.98) +
  scale_colour_manual(values = c("warn" = "#B2182B", "ok" = "#2166AC"),
                      guide = "none") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0.25, 9.6), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(6, 6, 6, 6))

cairo_pdf("figs/Fig1.pdf", width = 7.2, height = 6.6)
print(p)
dev.off()

png("figs/Fig1.png", width = 7.2, height = 6.6, units = "in", res = 400, type = "cairo")
print(p)
dev.off()

cat("Fig1 done\n")
