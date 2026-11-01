#!/usr/bin/env Rscript
# fig3_forest.R —— Fig 3：16 个 IgAN 风险等位 → 血液甲基化（随机效应）森林图
# 数据：results/fig3_forest.tsv（由 optionA_summary + optionA_direction 交叉得到）
# 规则：只纳入 GoDMC 随机效应 p<5e-8 的配对；固定效应显著但 RE 不显著者剔除
#      （CFH rs6677604/cg07209298：FE p≈0 但 RE p=0.13、I²=97.6，已剔除）

suppressPackageStartupMessages({ library(ggplot2); library(scales) })
ROOT <- Sys.getenv("IGAN_MR_ROOT", unset = ".")
setwd(ROOT)

d <- read.delim("results/fig3_forest.tsv", sep = "\t", stringsAsFactors = FALSE)
d$lbl <- sprintf("%s %s / %s", d$gene, d$rsid, d$cpg)
d$meta <- sprintf("I\u00b2=%s%% k=%s", round(as.numeric(d$het_isq)), d$n_studies)
band <- c("DEFA1 8p23.1" = "8p23", "HORMAD2 22q12.2" = "22q12",
          "ITGAM/ITGAX 16p11.2" = "16p11", "MADCAM1 19p13.3" = "19p13",
          "ST6GAL1 3q27.3" = "3q27", "TNFRSF13B 17p11.2" = "17p11",
          "CFH/CFHR 1q32" = "1q32")
d$band <- factor(band[d$region],
                 levels = band[unique(d$region[order(match(d$region,
                   c("DEFA1 8p23.1", "HORMAD2 22q12.2", "ITGAM/ITGAX 16p11.2",
                     "MADCAM1 19p13.3", "ST6GAL1 3q27.3", "TNFRSF13B 17p11.2",
                     "CFH/CFHR 1q32")))])])
d <- d[order(d$region, d$beta), ]
d$lbl <- factor(d$lbl, levels = rev(unique(d$lbl)))
d$lo <- d$beta - 1.96 * d$se
d$hi <- d$beta + 1.96 * d$se

theme_pub <- function(base = 9) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey90"),
          axis.text = element_text(colour = "grey15"),
          plot.title = element_text(face = "bold", size = base, hjust = 0),
          strip.background = element_rect(fill = "grey93", colour = NA),
          strip.text.y = element_text(face = "bold", size = base - 1.9, angle = 0))
}

p <- ggplot(d, aes(beta, lbl)) +
  geom_vline(xintercept = 0, linetype = "22", colour = "grey45", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.22, linewidth = 0.5,
                 colour = "grey35") +
  geom_point(aes(fill = beta > 0), shape = 21, size = 2.1, stroke = 0.45,
             colour = "grey15") +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2C6FA6"),
                    labels = c("risk allele lowers methylation",
                               "risk allele raises methylation"),
                    name = NULL) +
  geom_text(aes(x = 0.94, label = meta), hjust = 1, size = 2.05, colour = "grey40") +
  facet_grid(band ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(-0.85, 0.95), breaks = seq(-0.6, 0.6, 0.3)) +
  labs(x = "Change in blood methylation per copy of the IgAN risk allele (\u03b2, normalised units)",
       y = NULL,
       title = "IgAN risk alleles are methylation QTLs at five of seven susceptibility regions") +
  theme_pub() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        axis.text.y = element_text(size = 5.9),
        panel.spacing.y = unit(1.1, "mm"),
        plot.margin = margin(5, 30, 5, 5))

ggsave("figs/Fig3.pdf", p, width = 7.2, height = 4.9, device = cairo_pdf)
ggsave("figs/Fig3.png", p, width = 7.2, height = 4.9, dpi = 400)
cat(sprintf("Fig3 已写入 figs/：%d 个配对，%d 个区域\n", nrow(d), length(unique(d$region))))
