#!/usr/bin/env Rscript
# fig7_chromatin.R —— Fig 7：两个共定位 CpG 的染色质注释（hg38）
#
# 数据来源：
#   results/g4_ccre_annotation.tsv      ENCODE cCRE 区块（UCSC encodeCcreCombined）
#   results/g4_ccre_density_baseline.tsv  目标窗口密度 + 同染色体随机窗口基线
#   results/fig7_atac_lineage.tsv       覆盖 CpG 的 ChIP-Atlas ATAC 峰（血液细胞谱系拆分）
#   results/g4_atac_summary.tsv         ATAC 峰总数 / 血液数汇总
#
# 四个面板：
#   a  DEFA1 8p23.1 区域 cCRE 轨道（110 kb，22 个 cCRE；CpG 不在任何 cCRE 内）
#   b  HORMAD2 22q12.2 区域 cCRE 轨道（40 kb，44 个 cCRE；CpG 在 dELS EH38E2157595 内）
#   c  cCRE 密度 vs 同染色体随机窗口基线（第 50 vs 第 100 百分位）
#   d  覆盖 CpG 的血液 ATAC 峰谱系构成（0 vs 1,752）
#
# 输出：figs/Fig7.pdf、figs/Fig7.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
})

ROOT <- Sys.getenv("IGAN_MR_ROOT", unset = ".")
setwd(ROOT)
dir.create("figs", showWarnings = FALSE)

rd <- function(f) read.delim(file.path(ROOT, f), sep = "\t", stringsAsFactors = FALSE,
                             check.names = FALSE)

theme_pub <- function(base = 9) {
  theme_bw(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
      axis.text = element_text(colour = "grey15"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = base, hjust = 0,
                                margin = margin(b = 3)),
      plot.tag = element_text(face = "bold", size = base + 2),
      legend.key.size = unit(3.2, "mm"),
      strip.background = element_rect(fill = "grey93", colour = NA),
      strip.text = element_text(face = "bold", size = base - 0.5)
    )
}

CLASS_COL <- c("dELS" = "#4A7BB7", "pELS" = "#8FB4D9", "PLS" = "#B2182B",
               "CTCF-only" = "grey65")
CLASS_Y   <- c("dELS" = 3, "pELS" = 2, "PLS" = 4, "CTCF-only" = 1)

# ------------------------------------------------------------------ 读数据
cc <- rd("results/g4_ccre_annotation.tsv")
bl <- rd("results/g4_ccre_density_baseline.tsv")
atac_lin <- rd("results/fig7_atac_lineage.tsv")
atac_sum <- rd("results/g4_atac_summary.tsv")

make_track <- function(df, cpg_pos, title, xlab, tag, highlight_id = NULL,
                       gene_marks = NULL, note = NULL) {
  df <- df[order(-df$start), ]
  df$y <- CLASS_Y[df$ccre_class]
  df$fill <- df$ccre_class
  p <- ggplot() +
    # 基因标注（顶部）
    {if (!is.null(gene_marks))
      geom_segment(data = gene_marks,
                   aes(x = pos, xend = pos, y = 4.62, yend = 4.90),
                   colour = "#2166AC", linewidth = 0.5)} +
    {if (!is.null(gene_marks))
      geom_text(data = gene_marks, aes(x = pos, y = 5.06, label = name, hjust = hj),
                size = 2.2, colour = "#2166AC", fontface = "bold.italic")} +
    # cCRE 区块
    geom_rect(data = df,
              aes(xmin = start / 1e6, xmax = end / 1e6,
                  ymin = y - 0.30, ymax = y + 0.30, fill = fill),
              colour = "white", linewidth = 0.15) +
    {if (!is.null(highlight_id))
      geom_rect(data = df[df$ccre_id == highlight_id, ],
                aes(xmin = start / 1e6, xmax = end / 1e6,
                    ymin = y - 0.38, ymax = y + 0.38),
                fill = NA, colour = "black", linewidth = 0.6)} +
    scale_fill_manual(values = CLASS_COL, name = "cCRE class") +
    scale_y_continuous(breaks = c(1, 2, 3, 4),
                       labels = c("CTCF-only", "pELS", "dELS", "PLS"),
                       limits = c(0.4, 5.3)) +
    # CpG 位置
    geom_vline(xintercept = cpg_pos / 1e6, colour = "black", linetype = "dashed",
               linewidth = 0.5) +
    annotate("point", x = cpg_pos / 1e6, y = 0.55, shape = 25, fill = "#B2182B",
             colour = "black", size = 1.9) +
    labs(title = title, x = xlab, y = NULL, tag = tag) +
    theme_pub() +
    theme(panel.grid.major = element_blank(),
          axis.text.y = element_text(size = 7),
          legend.position = "none",
          plot.margin = margin(2, 6, 2, 4))
  if (!is.null(note))
    p <- p + annotate("text", x = cpg_pos / 1e6, y = 4.62, label = note,
                      size = 2.3, hjust = 1.02, colour = "grey25")
  p
}

# ------------------------------------------------------------------ 面板 a：DEFA1
d8 <- cc[cc$region == "DEFA1_8p23", ]
genes8 <- data.frame(pos = c(6969113, 6980094, 6999198) / 1e6,
                     name = c("DEFA10P TSS", "DEFA1 TSS", "DEFA1B TSS"),
                     hj = c(1.05, -0.05, -0.05))
pa <- make_track(d8, 6970292,
                 title = "DEFA1 8p23.1: 22 cCREs / 110 kb (0.20 per kb) - CpG in NO cCRE",
                 xlab = "chr8 position (Mb, hg38)", tag = "a",
                 gene_marks = genes8)

# ------------------------------------------------------------------ 面板 b：HORMAD2
d22 <- cc[cc$region == "HORMAD2_22q12", ]
genes22 <- data.frame(pos = 30207194 / 1e6, name = "LOC105372988 TSS", hj = -0.05)
pb <- make_track(d22, 30196307,
                 title = "HORMAD2 22q12.2: 44 cCREs / 40 kb (1.10 per kb) - CpG inside dELS",
                 xlab = "chr22 position (Mb, hg38)", tag = "b",
                 highlight_id = "EH38E2157595",
                 gene_marks = genes22,
                 note = "EH38E2157595 (dELS, score 587)")

# ------------------------------------------------------------------ 面板 c：密度基线
bl$group <- ifelse(bl$chrom == "chr8", "DEFA1 8p23.1", "HORMAD2 22q12.2")
bl$group <- factor(bl$group, levels = c("DEFA1 8p23.1", "HORMAD2 22q12.2"))
bl$type <- ifelse(bl$window_type == "TARGET", "Target window", "Random windows (same chr)")
bl$pct <- c("50th percentile", "100th percentile")[match(bl$chrom, c("chr8", "chr22"))]

pc <- ggplot(bl, aes(group, ccre_per_kb)) +
  geom_point(data = bl[bl$window_type == "RANDOM", ],
             aes(colour = type),
             position = position_jitter(width = 0.13, seed = 7),
             size = 1.5, shape = 16, alpha = 0.75) +
  geom_point(data = bl[bl$window_type == "TARGET", ],
             aes(colour = type), shape = 18, size = 3.4) +
  geom_text(data = bl[bl$window_type == "TARGET" & bl$chrom == "chr8", ],
            aes(label = sprintf("%.2f per kb\n(%s)", ccre_per_kb, pct)),
            nudge_x = 0.24, size = 2.3, colour = "#B2182B", hjust = 0) +
  geom_text(data = bl[bl$window_type == "TARGET" & bl$chrom == "chr22", ],
            aes(label = sprintf("%.2f per kb\n(%s)", ccre_per_kb, pct)),
            nudge_x = -0.24, nudge_y = 0.075, size = 2.3, colour = "#B2182B", hjust = 1) +
  scale_colour_manual(values = c("Target window" = "#B2182B",
                                 "Random windows (same chr)" = "grey55"),
                      name = NULL) +
  scale_y_continuous(limits = c(0, 1.32), breaks = seq(0, 1.2, 0.3)) +
  labs(title = "Regional cCRE density vs random baseline",
       x = NULL, y = "cCRE count per kb", tag = "c") +
  theme_pub() +
  theme(legend.position = "top", legend.text = element_text(size = 7))

# ------------------------------------------------------------------ 面板 d：血液 ATAC 谱系
atac_lin$lineage <- factor(atac_lin$lineage,
                           levels = c("T", "Myeloid", "B", "HSC/Progenitor", "other"))
atac_lin$lab <- ifelse(atac_lin$cpg == "cg00705004", "cg00705004 (DEFA1)",
                       "cg00539347 (HORMAD2)")
atac_lin$lab <- factor(atac_lin$lab, levels = c("cg00705004 (DEFA1)",
                                                "cg00539347 (HORMAD2)"))
LIN_COL <- c("T" = "#B2182B", "Myeloid" = "#C76A2A", "B" = "#2166AC",
             "HSC/Progenitor" = "#6A51A3", "other" = "grey75")

pd <- ggplot(atac_lin, aes(n, lab, fill = lineage)) +
  geom_col(width = 0.52, colour = "white", linewidth = 0.2) +
  annotate("text", x = 0, y = 2.47, label = "1,752 blood / 8,397 total",
           size = 2.5, hjust = 0, colour = "grey15") +
  annotate("text", x = 60, y = 1, label = "0 blood (3 total: H9 ESC only)",
           size = 2.5, hjust = 0, colour = "grey15") +
  scale_fill_manual(values = LIN_COL, name = "Blood lineage") +
  scale_x_continuous(limits = c(0, 2780), breaks = seq(0, 2000, 500),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "ATAC peaks covering the CpG (ChIP-Atlas, hg38)",
       x = "Number of ATAC peaks", y = NULL, tag = "d") +
  theme_pub() +
  theme(legend.position = "right", legend.text = element_text(size = 7))

# ------------------------------------------------------------------ 拼图输出
fig <- (pa / pb / (pc | pd)) + plot_layout(heights = c(0.95, 0.95, 1.35))

cairo_pdf("figs/Fig7.pdf", width = 7.2, height = 7.6)
print(fig)
dev.off()

png("figs/Fig7.png", width = 7.2, height = 7.6, units = "in", res = 400, type = "cairo")
print(fig)
dev.off()

cat("Fig7 done\n")
