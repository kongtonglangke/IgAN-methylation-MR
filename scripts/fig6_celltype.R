#!/usr/bin/env Rscript
# fig6_celltype.R —— Fig 6：22 个 CpG 在分选血细胞中的甲基化（GSE35069）
#
# 数据来源：
#   results/celltype_specificity.tsv   22 CpG × (7 分选细胞 + PBMC + 全血) 的平均甲基化 β
#
# 两个面板：
#   a  热图：22 CpG × 9 个细胞组分（按跨细胞类型极差降序排列）
#   b  极差条形图：每 CpG 的跨分选细胞极差，0.20 阈值虚线
#
# 输出：figs/Fig6.pdf、figs/Fig6.png

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

# ------------------------------------------------------------------ 数据
d <- rd("results/celltype_specificity.tsv")
# 排序：跨 7 个分选细胞（不含 PBMC/全血）的极差降序
sorted_cpg <- d$cpg[order(-d$range)]
d$cpg <- factor(d$cpg, levels = rev(sorted_cpg))

CELLS <- c("Neutrophil", "Granulocyte", "Eosinophil", "Monocyte",
           "B_cell", "CD4_T", "CD8_T", "NK", "PBMC", "WholeBlood")
CELL_LAB <- c("Neutro", "Gran", "Eosino", "Mono", "B cell", "CD4 T",
              "CD8 T", "NK", "PBMC", "Whole\nblood")

# 手工 melt（base R）
grid_cells <- expand.grid(idx = seq_len(nrow(d)), cell = CELLS, stringsAsFactors = FALSE)
grid_cells$beta <- mapply(function(i, cc) d[[cc]][i], grid_cells$idx, grid_cells$cell)
grid_cells$cpg <- factor(d$cpg[grid_cells$idx], levels = rev(sorted_cpg))
m <- grid_cells[, c("cpg", "cell", "beta")]
m$cell <- factor(m$cell, levels = CELLS, labels = CELL_LAB)

# 关键 CpG 行标签着色：共定位阳性 2 个（红/蓝）+ B 细胞特异 1 个（橙）
key <- c("cg00539347", "cg00705004", "cg03100639")
row_cols <- ifelse(sorted_cpg == "cg00539347", "#B2182B",
            ifelse(sorted_cpg == "cg00705004", "#2166AC",
            ifelse(sorted_cpg == "cg03100639", "#C76A2A", "grey20")))
row_face <- ifelse(sorted_cpg %in% key, "bold", "plain")

pa <- ggplot(m, aes(cell, cpg, fill = beta)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F6F2", high = "#B2182B",
                       midpoint = 0.5, limits = c(0, 1),
                       name = "Methylation\n(beta)") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Mean methylation per sorted blood cell type (GSE35069)",
       x = NULL, y = NULL, tag = "a") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6.2, colour = rev(row_cols),
                                   face = rev(row_face)),
        axis.text.x = element_text(size = 6.8, angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey40"),
        legend.position = "right")

# ------------------------------------------------------------------ 面板 b：极差
rg <- data.frame(cpg = sorted_cpg, range = d$range[order(-d$range)])
rg$cpg <- factor(rg$cpg, levels = rev(sorted_cpg))
rg$grp <- ifelse(rg$cpg == "cg00539347", "cg00539347 (coloc +)",
          ifelse(rg$cpg == "cg00705004", "cg00705004 (coloc +)",
          ifelse(rg$cpg == "cg03100639", "cg03100639 (B-cell specific)",
                 "other")))

pb <- ggplot(rg, aes(range, cpg, fill = grp)) +
  geom_col(width = 0.68, colour = "grey30", linewidth = 0.15) +
  geom_vline(xintercept = 0.20, linetype = "dashed", colour = "grey25",
             linewidth = 0.4) +
  annotate("text", x = 0.235, y = 1.6, label = "range = 0.20", size = 2.3,
           hjust = 0, colour = "grey25") +
  scale_fill_manual(values = c("cg00539347 (coloc +)" = "#B2182B",
                               "cg00705004 (coloc +)" = "#2166AC",
                               "cg03100639 (B-cell specific)" = "#C76A2A",
                               "other" = "grey75"),
                    guide = "none") +
  scale_x_continuous(limits = c(0, 0.72), breaks = seq(0, 0.6, 0.2)) +
  labs(title = "Cell-type range (max - min,\nsorted cells)",
       x = "Range of beta", y = NULL, tag = "b") +
  theme_pub() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 8))

# ------------------------------------------------------------------ 输出
fig <- pa + pb + plot_layout(widths = c(2.5, 1))

cairo_pdf("figs/Fig6.pdf", width = 7.2, height = 6.8)
print(fig)
dev.off()

png("figs/Fig6.png", width = 7.2, height = 6.8, units = "in", res = 400, type = "cairo")
print(fig)
dev.off()

cat("Fig6 done\n")
