#!/usr/bin/env Rscript
# fig8_gdia1.R —— Fig 8：O-糖基化通路（C1GALT1 / C1GALT1C1）的甲基化 QTL
#
# 数据来源：
#   results/gdia1_pathway_mqtl.tsv   ARIES(mQTLdb) 的 3 SNP × 5 时点 mQTL + Kiryluk 2017 Gd-IgA1 效应
#   results/celltype_c1galt1.tsv     GSE35069，7 个 CpG × 7 类分选血细胞 + PBMC + 全血
#
# 四个面板：
#   a  rs13226913 在 5 个发育时点对 cg04337386 的效应（森林图）
#   b  三个 SNP 统一到 Gd-IgA1 升高等位后的甲基化效应方向
#   c  细胞类型甲基化分布（cg04337386 近乎完全未甲基化）
#   d  方向链示意：mQTL / eQTL / Gd-IgA1 的同向关系不支持甲基化中介
#
# 输出：figs/Fig8.pdf、figs/Fig8.png
#
# 注意：mQTL 的 SE 未在公开记录中给出，由 beta 与双侧 p 反算（z = qnorm(p/2, lower=FALSE)），
#       仅用于 95% CI 的图形展示，不用于任何推断。beta 为 M-value 尺度。

suppressPackageStartupMessages({
  library(ggplot2)
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

COL_UP   <- "#B2182B"   # 甲基化升高
COL_DOWN <- "#2C6FA6"   # 甲基化降低

# ------------------------------------------------------------------ 读数据
d <- rd("results/gdia1_pathway_mqtl.tsv")

# 统一到"Gd-IgA1 升高等位"：inc 等位 == a1 时保留原 beta，否则取相反数
d$beta_inc <- ifelse(d$gdia1_inc_allele == d$a1, d$mqtl_beta_a1, -d$mqtl_beta_a1)
# SE 由双侧 p 反算（仅作图用）
d$zstat <- qnorm(d$mqtl_p / 2, lower.tail = FALSE)
d$se    <- abs(d$mqtl_beta_a1) / d$zstat
d$lo    <- d$beta_inc - 1.96 * d$se
d$hi    <- d$beta_inc + 1.96 * d$se

tp_levels <- c("Birth", "Childhood", "Adolescence", "MiddleAge", "Pregnancy")
tp_labels <- c("Birth", "Childhood", "Adolescence", "Middle age", "Pregnancy")

# ------------------------------------------------------------------ Panel a
da <- d[d$snp == "rs13226913", ]
da$tp <- factor(da$timepoint, levels = tp_levels, labels = tp_labels)
da <- da[order(da$tp), ]
da$plab <- sprintf("%.1e", da$mqtl_p)

pa <- ggplot(da, aes(beta_inc, tp)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "grey55", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.55,
                colour = COL_DOWN, orientation = "y") +
  geom_point(size = 2.1, colour = COL_DOWN, shape = 21, fill = COL_DOWN) +
  geom_text(aes(x = 0.040, label = plab), hjust = 0, size = 2.3, colour = "grey35") +
  scale_x_continuous(limits = c(-0.50, 0.098), breaks = seq(-0.5, 0, 0.1)) +
  labs(x = expression(paste("Methylation change per T allele (", beta, ", normalised units)")),
       y = NULL,
       title = "rs13226913 \u2192 cg04337386, five ages",
       tag = "a") +
  theme_pub() +
  theme(panel.grid.major.y = element_blank())

# ------------------------------------------------------------------ Panel b
# 每个 SNP 取 p 最小的时点作为代表
rep_rows <- do.call(rbind, lapply(split(d, d$snp), function(x) x[which.min(x$mqtl_p), ]))
rep_rows$locus <- c(rs1008897 = "C1GALT1", rs13226913 = "C1GALT1",
                    rs5910940 = "C1GALT1C1")[rep_rows$snp]
rep_rows$ord <- rep_rows$beta_inc
rep_rows <- rep_rows[order(rep_rows$ord), ]
rep_rows$lab <- factor(sprintf("%s\n%s  \u00b7  %s", rep_rows$snp, rep_rows$locus, rep_rows$cpg),
                       levels = sprintf("%s\n%s  \u00b7  %s", rep_rows$snp, rep_rows$locus, rep_rows$cpg))
rep_rows$dir <- ifelse(rep_rows$beta_inc < 0, "demethylating", "methylating")

pb <- ggplot(rep_rows, aes(beta_inc, lab)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "grey55", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.55,
                colour = "grey35", orientation = "y") +
  geom_point(aes(fill = dir), size = 2.9, shape = 21, colour = "grey20", stroke = 0.35) +
  geom_text(aes(x = hi + 0.035, label = sprintf("%.2f", beta_inc)),
            hjust = 0, size = 2.3, colour = "grey20") +
  scale_fill_manual(values = c(demethylating = COL_DOWN, methylating = COL_UP), guide = "none") +
  scale_x_continuous(limits = c(-0.50, 0.26), breaks = seq(-0.4, 0.2, 0.2)) +
  labs(x = expression(paste("Methylation change per Gd-IgA1-raising allele (", beta, ")")),
       y = NULL,
       title = "Aligned effect direction",
       tag = "b") +
  theme_pub() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 6.2, lineheight = 0.95),
        axis.title.x = element_text(size = 7.6))

# ------------------------------------------------------------------ Panel c
ct <- rd("results/celltype_c1galt1.tsv")
cell_cols <- c("Granulocyte", "Eosinophil", "Monocyte", "B_cell",
               "CD4_T", "CD8_T", "NK", "PBMC", "WholeBlood")
ct_long <- do.call(rbind, lapply(cell_cols, function(cc) {
  data.frame(cpg = ct$cpg, cell = cc, beta = ct[[cc]],
             stringsAsFactors = FALSE)
}))
ct_long$cell <- factor(ct_long$cell, levels = cell_cols,
                       labels = c("Granulocyte", "Eosinophil", "Monocyte", "B cell",
                                  "CD4 T", "CD8 T", "NK", "PBMC", "Whole blood"))
# CpG 按甲基化水平排序，cg04337386 自然落到最低端
ord <- ct[order(ct$Granulocyte + ct$B_cell), ]
ct_long$cpg <- factor(ct_long$cpg, levels = ord$cpg)
cpg_lab <- setNames(sprintf("%s   (range %.3f)", ct$cpg, ct$range), ct$cpg)
# 9 列格宽不足以容纳三位小数；数值改由图注给出，格内只靠颜色表达
pc <- ggplot(ct_long, aes(cell, cpg, fill = beta)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  scale_fill_gradient2(low = "#2C6FA6", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0.4, limits = c(0, 1),
                       name = "Methylation\n(beta value)") +
  scale_y_discrete(labels = cpg_lab) +
  labs(x = NULL, y = NULL,
       title = "C1GALT1 locus in sorted blood cells",
       tag = "c") +
  theme_pub() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 1, size = 6.2),
        axis.text.y = element_text(size = 6.1),
        legend.key.height = unit(7, "mm"),
        legend.title = element_text(size = 6.4),
        legend.text = element_text(size = 6.2))

# ------------------------------------------------------------------ Panel d
dl <- data.frame(
  x = c(1.0, 5.0, 5.0, 9.2),
  y = c(5.0, 8.2, 1.8, 5.0),
  lab = c("rs13226913\nT allele",
          "cg04337386\nmethylation \u2193\n(\u03b2 = \u22120.38)",
          "C1GALT1 mRNA \u2193\n(cis-eQTL\np = 3.9e-23)",
          "Serum Gd-IgA1 \u2191\n(+0.22 SD,\np = 3.2e-11)"),
  stringsAsFactors = FALSE
)

pd <- ggplot() +
  annotate("segment", x = 1.9, y = 5.6, xend = 4.1, yend = 7.7,
           linewidth = 0.55, colour = COL_DOWN,
           arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  annotate("segment", x = 1.9, y = 4.4, xend = 4.1, yend = 2.3,
           linewidth = 0.55, colour = COL_DOWN,
           arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  annotate("segment", x = 5.0, y = 7.0, xend = 5.0, yend = 2.6,
           linewidth = 0.5, colour = "grey45", linetype = "22") +
  annotate("segment", x = 5.9, y = 7.7, xend = 8.1, yend = 5.9,
           linewidth = 0.55, colour = "grey25",
           arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  annotate("segment", x = 5.9, y = 2.3, xend = 8.1, yend = 4.1,
           linewidth = 0.55, colour = "grey25",
           arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  annotate("label", x = 5.35, y = 5.0, label = "same direction\n(no mediation)",
           size = 2.0, colour = "grey35", hjust = 0, lineheight = 0.95,
           label.size = 0, fill = "white", label.padding = unit(1.0, "mm")) +
  annotate("text", x = 3.0, y = 7.15, label = "mQTL", size = 2.2, colour = COL_DOWN) +
  annotate("text", x = 3.0, y = 2.85, label = "eQTL", size = 2.2, colour = COL_DOWN) +
  geom_label(data = dl, aes(x, y, label = lab), size = 2.1, lineheight = 0.92,
             label.size = 0.18, label.padding = unit(1.3, "mm"),
             fill = c("grey96", "#EAF1F8", "#EAF1F8", "#FBEDEC"),
             colour = "grey15") +
  coord_cartesian(xlim = c(0.2, 10.8), ylim = c(0.4, 10.0)) +
  labs(x = NULL, y = NULL,
       title = "Methylation and expression move together",
       tag = "d") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0,
                                  margin = margin(b = 3)),
        plot.tag = element_text(face = "bold", size = 11),
        plot.margin = margin(4, 6, 4, 4))

# ------------------------------------------------------------------ 输出
lay <- grid.layout(2, 2,
                   widths  = unit(c(1, 1), "null"),
                   heights = unit(c(1, 1.16), "null"))

cairo_pdf("figs/Fig8.pdf", width = 7.2, height = 6.4)
grid.newpage()
pushViewport(viewport(layout = lay))
print(pa, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pb, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pc, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(pd, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

png("figs/Fig8.png", width = 7.2, height = 6.4, units = "in", res = 400, type = "cairo")
grid.newpage()
pushViewport(viewport(layout = lay))
print(pa, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pb, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pc, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(pd, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

cat("Fig8 已写入 figs/\n")
cat(sprintf("a: %d 时点，beta_inc %.3f ~ %.3f\n", nrow(da), min(da$beta_inc), max(da$beta_inc)))
cat(sprintf("b: %s\n", paste(rep_rows$snp, rep_rows$timepoint, sprintf("%.2f", rep_rows$beta_inc),
                              sep = "@", collapse = ", ")))
cat(sprintf("c: %d CpG x %d 类型；cg04337386 range = %.4f\n",
            nrow(ct), length(cell_cols), ct$range[ct$cpg == "cg04337386"]))
