#!/usr/bin/env Rscript
# fig5_defa1_cnv.R —— Fig 5：DEFA1 位点的拷贝数 / 重复 / 探针特异性检查（hg19）
#
# 数据来源：
#   results/fig5_segdup.tsv          UCSC segmental duplication（chr8 6.66-6.97 Mb, >=2 kb & >=90%）
#   results/cnv_confound_ref.json    4 个 CpG 的 segdup 归属 + 50-nt k-mer 唯一性
#   results/defa1_cnv_signature.json 剂量签名（lead SNP 显著的 CpG 数 / 各自 lead SNP 数）
#   results/coloc_cnv_sensitivity.tsv  去重复块 / 去 driver SNP 后的 coloc 敏感性
#   data/gene_coords_grch37.tsv      DEFA1 基因坐标（hg19）
#
# 三个面板：
#   a  位点图：片段重复块（y=比对一致性）、19 kb 串联重复、DUP 结构变异、4 个 CpG、DEFA1 基因
#   b  剂量签名：lead SNP 显著 CpG 占比 vs 各 CpG 独立 lead SNP 占比（三基因对照）
#   c  coloc 敏感性：去掉重复块 / driver SNP ±2 kb 后 PP.H4（cg00539347、cg13682187 为对照）
#
# 输出：figs/Fig5.pdf、figs/Fig5.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
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
sd  <- rd("results/fig5_segdup.tsv")
sd$start <- sd$start / 1e6
sd$end   <- sd$end / 1e6

sig <- jsonlite::fromJSON(file.path(ROOT, "results/defa1_cnv_signature.json"))
cnv <- jsonlite::fromJSON(file.path(ROOT, "results/cnv_confound_ref.json"))
sens <- rd("results/coloc_cnv_sensitivity.tsv")

# 关键位点（hg19）
REPEAT <- c(6.839960, 6.878169)          # 19 kb DEFA1/DEFA1B/DEFA3 串联重复，99.52%
DUP_SV <- c(6.879829, 6.920559)          # 1000G phase3 唯一 SV（DUP，EUR_AF=0）
DEFA1  <- c(6.835171, 6.837602)          # DEFA1 基因（hg19，负链）
CPG <- data.frame(
  cpg = c("cg00705004", "cg07297068", "cg19304352", "cg05969298"),
  pos = c(6827814, 6691691, 6796176, 6924737) / 1e6,
  inside = c(FALSE, FALSE, FALSE, TRUE)
)
CPG$label <- paste0("italic('", CPG$cpg, "')")

# ------------------------------------------------------------------ 面板 a：位点图
p_region <- data.frame(
  x = c(REPEAT[1], DUP_SV[1]),
  xend = c(REPEAT[2], DUP_SV[2]),
  y = c(1.006, 0.902), yend = c(1.006, 0.902),
  what = c("Tandem repeat 19 kb (99.52%)", "DUP SV (1000G, EUR AF = 0)")
)

pa <- ggplot() +
  # 重复块：水平线段，y=一致性
  geom_segment(data = sd, aes(x = start, xend = end, y = identity, yend = identity,
                              colour = identity), linewidth = 1.6, alpha = 0.85) +
  scale_colour_gradientn(colours = c("#F4E9D4", "#E8B96B", "#C76A2A"),
                         limits = c(0.90, 1.0), name = "Duplication identity",
                         guide = guide_colorbar(breaks = c(0.90, 0.95, 1.00),
                                                barwidth = unit(3.6, "cm"),
                                                barheight = unit(0.45, "cm"))) +
  # 19 kb 串联重复与 DUP SV
  geom_segment(data = p_region,
               aes(x = x, xend = xend, y = y, yend = yend, linetype = what),
               linewidth = 1.1, colour = c("#8C510A", "#5F5F5F"),
               show.legend = TRUE) +
  scale_linetype_manual(values = c("longdash", "solid"), guide = "none") +
  # DEFA1 基因
  annotate("rect", xmin = DEFA1[1], xmax = DEFA1[2], ymin = 0.896, ymax = 0.900,
           fill = "#2166AC") +
  annotate("text", x = (DEFA1[1] + DEFA1[2]) / 2, y = 0.890, label = "DEFA1",
           size = 2.2, colour = "#2166AC", fontface = "bold") +
  # CpG 位置
  geom_vline(data = CPG, aes(xintercept = pos),
             colour = ifelse(CPG$cpg == "cg00705004", "#B2182B", "grey55"),
             linetype = ifelse(CPG$cpg == "cg00705004", "solid", "dashed"),
             linewidth = ifelse(CPG$cpg == "cg00705004", 0.55, 0.35)) +
  annotate("point", x = CPG$pos, y = 1.006,
           shape = c(18, 16, 16, 16), size = c(2.6, 1.8, 1.8, 1.8),
           colour = c("#B2182B", "grey35", "grey35", "#C76A2A")) +
  annotate("text", x = CPG$pos[1] - 0.004, y = 1.0145, label = "cg00705004",
           size = 2.3, colour = "#B2182B", fontface = "bold", hjust = 1) +
  annotate("text", x = CPG$pos[4] + 0.003, y = 1.0145, label = "cg05969298",
           size = 2.3, colour = "#C76A2A", hjust = 0) +
  annotate("text", x = CPG$pos[2], y = 0.883, label = "cg07297068", size = 2.2,
           colour = "grey35") +
  annotate("text", x = CPG$pos[3], y = 0.883, label = "cg19304352", size = 2.2,
           colour = "grey35") +
  annotate("text", x = 6.72, y = 1.0205, size = 2.4, hjust = 0, colour = "grey25",
           label = "Probe uniqueness: 50-nt k-mer, <=2 mismatches -> 0 paralogous copies for all 4 CpGs") +
  coord_cartesian(ylim = c(0.88, 1.022), xlim = c(6.66, 6.97), clip = "off") +
  scale_y_continuous(limits = c(0.88, 1.022), breaks = c(0.90, 0.95, 1.00)) +
  scale_x_continuous(breaks = seq(6.7, 6.95, 0.05)) +
  labs(title = "Defensin cluster, chr8p23.1 (hg19)",
       x = "Chromosome 8 position (Mb)", y = "Segmental duplication identity",
       tag = "a") +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.box.margin = margin(t = -2),
        plot.margin = margin(12, 8, 4, 4))

# ------------------------------------------------------------------ 面板 b：剂量签名
dose <- data.frame(
  gene = rep(c("DEFA1", "HORMAD2", "CFH"), 2),
  metric = rep(c("CpGs significant\nfor lead SNP", "CpGs with distinct\nown-lead SNP"), each = 3),
  pct = c(100 * sig$n_sig / sig$n_cpg, 100 * sig$n_distinct_lead / sig$n_cpg),
  lab = c(sprintf("%d/%d", sig$n_sig, sig$n_cpg),
          sprintf("%d/%d", sig$n_distinct_lead, sig$n_cpg))
)
dose$gene <- factor(dose$gene, levels = c("DEFA1", "HORMAD2", "CFH"))
dose$metric <- factor(dose$metric, levels = c("CpGs significant\nfor lead SNP",
                                              "CpGs with distinct\nown-lead SNP"))
# 手工置信区间（Wilson 粗略）：用点式条形即可，不画 CI（描述性）
pb <- ggplot(dose, aes(gene, pct, fill = metric)) +
  geom_col(position = position_dodge(0.72), width = 0.62, colour = "grey25",
           linewidth = 0.2) +
  geom_text(aes(label = lab), position = position_dodge(0.72), vjust = -0.45,
            size = 2.3, colour = "grey15") +
  scale_fill_manual(values = c("#C76A2A", "#4A7BB7"), name = NULL) +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "Dose-signature test: no single driver",
       x = NULL, y = "Share of CpGs in region", tag = "b") +
  theme_pub() +
  theme(legend.position = "top",
        legend.text = element_text(size = 7))

# ------------------------------------------------------------------ 面板 c：coloc 敏感性
sens$cpg <- factor(sens$cpg,
                   levels = c("cg00705004", "cg00539347", "cg13682187"),
                   labels = c("cg00705004", "cg00539347", "cg13682187"))
sens$snpset <- factor(sens$snpset,
                      levels = c("all", "outside_dup", "outside_dup_2k"),
                      labels = c("All SNPs", "- repeat block", "- repeat - driver\n±2 kb"))
pc <- ggplot(sens, aes(snpset, H4, fill = cpg, group = cpg)) +
  geom_col(position = position_dodge(0.78), width = 0.68, colour = "grey25",
           linewidth = 0.2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "#B2182B",
             linewidth = 0.4) +
  annotate("text", x = 3.56, y = 1.075, label = "PP.H4 > 0.8", size = 2.3,
           colour = "#B2182B", hjust = 1) +
  geom_text(aes(label = sprintf("%.3f", H4)),
            position = position_dodge(0.78), vjust = -0.4, size = 2.1,
            colour = "grey15") +
  scale_fill_manual(values = c("#B2182B", "#4A7BB7", "grey65"), name = NULL) +
  scale_y_continuous(limits = c(0, 1.14), breaks = seq(0, 1, 0.25)) +
  labs(title = "Colocalization robust to repeat removal",
       x = NULL, y = "PP.H4 (shared causal variant)", tag = "c") +
  theme_pub() +
  theme(legend.position = "top", legend.text = element_text(size = 7))

# ------------------------------------------------------------------ 拼图输出
fig <- pa / (pb | pc) + plot_layout(heights = c(1.15, 1))

cairo_pdf("figs/Fig5.pdf", width = 7.2, height = 6.6)
print(fig)
dev.off()

png("figs/Fig5.png", width = 7.2, height = 6.6, units = "in", res = 400, type = "cairo")
print(fig)
dev.off()

cat("Fig5 done\n")
