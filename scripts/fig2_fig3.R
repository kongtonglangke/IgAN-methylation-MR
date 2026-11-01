#!/usr/bin/env Rscript
# fig2_fig3.R —— 生成 Fig 2（恒等式散点 + Ⅰ类错误热图）与 Fig 3（coloc 后验堆叠图）
# 数据来源：
#   results/fig2a_identity.tsv      1699 个 CpG-SNP 配对的 z_MR 与 z_GWAS
#   results/sim_type1_table.tsv     固定种子 20261030 的 Ⅰ类错误模拟（若存在）
#   results/coloc_summary.tsv       19 个 CpG 的 coloc 后验
# 输出：figs/Fig2.pdf/.png、figs/Fig3.pdf/.png

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

# ---------------------------------------------------------------- Fig 2A
id <- rd("results/fig2a_identity.tsv")
dz <- abs(id$z_mr - id$z_gwas)
lab_txt <- sprintf("n = %s CpG\u2013SNP pairs\nmax |\u0394z| = %.1e\nmedian |\u0394z| = %.1e",
                   format(nrow(id), big.mark = ","), max(dz), median(dz))

pA <- ggplot(id, aes(z_gwas, z_mr)) +
  geom_abline(slope = 1, intercept = 0, colour = "#C0392B", linewidth = 0.45,
              linetype = "22") +
  geom_point(alpha = 0.18, size = 0.65, colour = "#1F3B73", stroke = 0) +
  annotate("label", x = -3.85, y = 3.6, label = lab_txt, hjust = 0, vjust = 1,
           size = 2.5, linewidth = 0.18, label.padding = unit(1.5, "mm"),
           fill = "white", colour = "grey20") +
  coord_equal(xlim = c(-4.2, 4.0), ylim = c(-4.2, 4.0)) +
  labs(x = expression(paste("GWAS ", italic(z), " for IgAN")),
       y = expression(paste("Wald-ratio ", italic(z), " of methylation MR")),
       title = "Wald ratio z equals GWAS z",
       tag = "a") +
  theme_pub()

# ---------------------------------------------------------------- Fig 2B
simf <- file.path(ROOT, "results/sim_type1_table.tsv")
if (file.exists(simf)) {
  sim <- rd("results/sim_type1_table.tsv")
  sim$method <- factor(sim$method,
                       levels = c("sumz", "pct", "fisher", "chisq", "minp", "pctcr"),
                       labels = c("Stouffer (\u03a3z/\u221ak)", "Naive z (control set)", "Fisher",
                                  "\u03c7\u00b2 (sum of z\u00b2)", "min-p (\u0160id\u00e1k)",
                                  "Gene-cluster t"))
  sim$rho_f <- factor(sprintf("\u03c1 = %.1f", sim$rho))
  sim$tile <- sprintf("%.3f", sim$err05)
  sim$k_f <- factor(ifelse(sim$k == 1699,
                           "Real CpG density (k = 1,699)",
                           "3\u00d7 scaled density (k = 5,097)"),
                    levels = c("Real CpG density (k = 1,699)",
                               "3\u00d7 scaled density (k = 5,097)"))
  pB <- ggplot(sim, aes(rho_f, method)) +
    geom_tile(aes(fill = err05), colour = "white", linewidth = 0.5) +
    geom_text(aes(label = tile, colour = err05 > 0.45), size = 2.1) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey10"), guide = "none") +
    scale_fill_gradient2(low = "#2C6FA6", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0.05, limits = c(0, 0.9),
                         name = "Type I error\n(\u03b1 = 0.05)") +
    facet_wrap(~ k_f, ncol = 1) +
    labs(x = "Within-gene correlation between instruments",
         y = NULL,
         title = "Type I error (2,000 replicates)") +
    theme_pub() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 6.6, lineheight = 0.88),
          axis.text.x = element_text(size = 6.8),
          strip.text = element_text(size = 6.8, face = "bold"),
          strip.background = element_rect(fill = "grey92", colour = NA),
          legend.position = "right",
          legend.key.height = unit(9, "mm"))
} else {
  pB <- ggplot() +
    annotate("text", x = 0, y = 0, label = "type-I simulation table not found",
             size = 3, colour = "grey40") +
    theme_void() + labs(tag = "b")
}

pdf("figs/Fig2.pdf", width = 7.2, height = 4.15, useDingbats = FALSE)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(1, 1.06), "null"))))
print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

png("figs/Fig2.png", width = 7.2, height = 4.15, units = "in", res = 400)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(1, 1.06), "null"))))
print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

# ---------------------------------------------------------------- Fig 3
co <- rd("results/coloc_summary.tsv")
co <- co[co$status == "ok", ]
reg_lab <- c(
  DEFA1_8p23       = "DEFA1  8p23.1",
  HORMAD2_22q12    = "HORMAD2/LOC105372988  22q12.2",
  CFH_CFHR_1q32    = "CFH/CFHR  1q32.2",
  ITGAM_ITGAX_16p11 = "ITGAM/ITGAX  16p11.2",
  MADCAM1_19p13    = "MADCAM1  19p13.3",
  ST6GAL1_3q27     = "ST6GAL1  3q27.3",
  TNFRSF13B_17p11  = "TNFRSF13B  17p11.2"
)
co$region_lab <- factor(reg_lab[co$region], levels = reg_lab[names(reg_lab) %in% co$region])
co <- co[order(co$region_lab, -co$PP.H4), ]
co$cpg <- factor(co$cpg, levels = rev(unique(co$cpg)))

long <- do.call(rbind, lapply(c("H0", "H1", "H2", "H3", "H4"), function(h) {
  data.frame(cpg = co$cpg, region_lab = co$region_lab,
             hyp = h, pp = co[[paste0("PP.", h)]])
}))
long$hyp <- factor(long$hyp, levels = c("H0", "H1", "H2", "H3", "H4"))

pal <- c(H0 = "grey88", H1 = "#9ECAE1", H2 = "#6BAED6", H3 = "#F4A582", H4 = "#B2182B")

p3 <- ggplot(long, aes(pp, cpg, fill = hyp)) +
  geom_col(width = 0.72, colour = NA) +
  geom_text(data = co[co$PP.H4 > 0.5, ],
            aes(x = PP.H4 - 0.015, y = cpg, label = sprintf("H4 = %.3f", PP.H4)),
            inherit.aes = FALSE, hjust = 1, size = 2.35, colour = "white") +
  facet_grid(region_lab ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = pal,
                    labels = c("H0: none", "H1: mQTL only", "H2: IgAN only",
                               "H3: distinct", "H4: shared"),
                    name = NULL) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.001)) +
  labs(x = "Posterior probability", y = NULL,
       title = "Colocalisation of blood mQTL and IgAN risk at seven susceptibility regions") +
  theme_pub() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(4.2, "mm"),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 6.4),
        strip.text.y = element_text(angle = 0, hjust = 0, size = 6.4),
        panel.spacing.y = unit(1.0, "mm"))

ggsave("figs/Fig4.pdf", p3, width = 6.6, height = 5.0, device = cairo_pdf)
ggsave("figs/Fig4.png", p3, width = 6.6, height = 5.0, dpi = 400)

cat("Fig2/Fig4 已写入 figs/\n")
cat(sprintf("Fig2A: n=%d, max|dz|=%.2e\n", nrow(id), max(dz)))
cat(sprintf("Fig4 : %d 个 CpG, %d 个区域\n", nrow(co), length(unique(co$region))))
