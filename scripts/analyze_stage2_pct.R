# Stage 2b: the calibrated pathway-level statistic applied to the real data.
#
# The type-I simulation (scripts/sim_type1.py) showed that every independence
# based combiner (Fisher, chi-square on mean z^2, sum z/sqrt(k), and even the
# naive control-set z) rejects 27-88 % of the time at alpha=0.05 once the
# within-gene instrument correlation exceeds 0.3, while the gene-clustered
# t-test stays at 4.6-8.6 %. Only the clustered statistic is used here.
suppressMessages(library(data.table))
setwd("F:/longxiawenjian/qiguanzhou/igan_epi_mr")

d <- fread("results/mr_stage1b_cpg_igan.tsv", sep = "\t", data.table = FALSE)
d <- d[d$status == "ok", ]
d <- d[!duplicated(d$cpg), ]
d$p <- as.numeric(d$pval)
d$b <- as.numeric(d$b)
d$z <- qnorm(d$p / 2, lower.tail = FALSE) * sign(d$b)   # z along the MR effect
d$zabs <- abs(qnorm(d$p / 2, lower.tail = FALSE))

cat("=== pathway members ===\n")
cat(sprintf("CpGs %d | genes %d | instruments %d\n",
            nrow(d), length(unique(d$gene)), length(unique(d$rsid))))

# gene-level clusters
gm <- do.call(rbind, lapply(split(d, d$gene), function(x) {
  data.frame(gene = x$gene[1], n = nrow(x),
             mz = mean(x$zabs), sdz = sd(x$zabs))
}))
gm <- gm[order(-gm$n), ]
print(format(gm, digits = 3), row.names = FALSE)

G <- nrow(gm)
# null expectation of |z| under no effect: E|N(0,1)| = sqrt(2/pi)
null_absz <- sqrt(2 / pi)
mg <- mean(gm$mz)
sg <- sd(gm$mz)
se <- sg / sqrt(G)
t <- (mg - null_absz) / se
p_clust <- 2 * pt(-abs(t), df = G - 1)
cat(sprintf("\n=== gene-clustered pathway test ===\n"))
cat(sprintf("null E|z| = %.3f\n", null_absz))
cat(sprintf("observed mean |z| across genes %.3f (sd %.3f, G=%d)\n", mg, sg, G))
cat(sprintf("cluster-robust SE %.3f -> t=%.2f, df=%d, p=%.2f\n", se, t, G - 1, p_clust))
cat(sprintf("95%% CI for the mean |z|: %.3f - %.3f\n",
            mg - qt(.975, G - 1) * se, mg + qt(.975, G - 1) * se))

# mean z^2 (two-sided polygenic signal), null = 1
gm2 <- do.call(rbind, lapply(split(d, d$gene), function(x)
  data.frame(gene = x$gene[1], nz2 = mean(x$zabs^2))))
t2 <- (mean(gm2$nz2) - 1) / (sd(gm2$nz2) / sqrt(nrow(gm2)))
cat(sprintf("mean z^2 %.3f (null 1) -> t=%.2f, p=%.2f\n",
            mean(gm2$nz2), t2, 2 * pt(-abs(t2), nrow(gm2) - 1)))

# What pathway-average effect could this design have detected?
cat("\n=== sensitivity of the pathway test ===\n")
for (pw in c(0.8, 0.9)) {
  need <- (qt(1 - 0.025, G - 1) + qt(pw, G - 1)) * se
  cat(sprintf("  %d%% power: mean |z| must exceed %.3f (i.e. +%.0f%% relative)\n",
              pw * 100, null_absz + need, 100 * need / null_absz))
}
# a uniform per-SD log-OR beta shifts every z by delta = beta/SE_MR; invert the
# folded-normal mean to translate a needed |z| shift into that beta
eabs <- function(delta) 2 * dnorm(delta) + delta * (1 - 2 * pnorm(-delta))
need <- (qt(1 - 0.025, G - 1) + qt(0.8, G - 1)) * se
delta <- uniroot(function(d) eabs(d) - null_absz - need, c(0, 5))$root
cat(sprintf("  that +%.3f in mean |z| needs a uniform z shift of %.2f,\n", need, delta))
cat(sprintf("  i.e. a per-SD log-OR of %.3f (OR %.2f) at every CpG\n",
            delta * median(as.numeric(d$se)), exp(delta * median(as.numeric(d$se)))))

cat("\n=== how many independent signals are there really? ===\n")
cat(sprintf("  CpGs %d -> distinct lead SNPs %d (%.1f CpGs per SNP)\n",
            nrow(d), length(unique(d$rsid)), nrow(d) / length(unique(d$rsid))))
cat(sprintf("  distinct lead SNPs per gene: median %.0f, range %d-%d\n",
            median(gm$n), min(gm$n), max(gm$n)))
