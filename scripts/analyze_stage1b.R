# Stage 1b summary: CpG methylation -> IgA nephropathy (Sakaue 2021, GCST90018866)
# Pure summary-statistics work, so R is used for the distributional tests.
suppressMessages(library(data.table))

setwd("F:/longxiawenjian/qiguanzhou/igan_epi_mr")
d <- fread("results/mr_stage1b_cpg_igan.tsv", sep = "\t", data.table = FALSE)
cat("rows:", nrow(d), "\n")
print(table(d$status))

ok <- d[d$status == "ok", ]
ok <- ok[!duplicated(ok$cpg), ]          # a CpG can map to overlapping genes
ok$p <- as.numeric(ok$pval)
ok$b <- as.numeric(ok$b)
ok$se <- as.numeric(ok$se)
ok$bx <- as.numeric(ok$bx)
ok$pout <- as.numeric(ok$p_out)
ok$F <- as.numeric(ok$lead_F)

cat("\n=== cohort of analysed CpGs ===\n")
cat(sprintf("unique CpGs %d | genes %d | independent lead SNPs %d\n",
            nrow(ok), length(unique(ok$gene)), length(unique(ok$rsid))))
cat(sprintf("median instrument F %.0f (IQR %.0f-%.0f)\n",
            median(ok$F), quantile(ok$F, .25), quantile(ok$F, .75)))

# the Wald ratio with one instrument has z identical to the SNP-outcome z,
# so the MR p-value is the GWAS p-value rescaled; state it explicitly
cat(sprintf("\nmax |z_MR - z_outcome| = %.3g  (expected 0 by construction)\n",
            max(abs(qnorm(ok$p / 2) - qnorm(ok$pout / 2)), na.rm = TRUE)))

cat("\n=== tail calibration (unique CpGs) ===\n")
n <- nrow(ok)
for (th in c(0.05, 0.01, 1e-3, 1e-4)) {
  k <- sum(ok$p < th)
  cat(sprintf("  p<%-8g  observed %4d  expected %7.1f  ratio %5.2f  binom.p %.3f\n",
              th, k, th * n, k / (th * n), binom.test(k, n, th)$p.value))
}

# collapse to one row per independent instrument (strongest CpG per lead SNP)
ok <- ok[order(ok$p), ]
ins <- ok[!duplicated(ok$rsid), ]
cat("\n=== tail calibration (independent instruments) ===\n")
m <- nrow(ins)
for (th in c(0.05, 0.01, 1e-3, 1e-4)) {
  k <- sum(ins$p < th)
  cat(sprintf("  p<%-8g  observed %4d  expected %7.1f  ratio %5.2f  binom.p %.3f\n",
              th, k, th * m, k / (th * m), binom.test(k, m, th)$p.value))
}

cat("\n=== uniformity ===\n")
ks <- ks.test(ins$p, "punif")
cat(sprintf("  KS D=%.4f p=%.3f (independent instruments, n=%d)\n",
            ks$statistic, ks$p.value, m))
chi <- qchisq(ins$p, 1, lower.tail = FALSE)
cat(sprintf("  lambda_GC = %.3f\n", median(chi) / qchisq(0.5, 1)))

cat("\n=== multiplicity ===\n")
cat(sprintf("  Bonferroni alpha/n = %.2e ; CpGs passing: %d\n",
            0.05 / n, sum(ok$p < 0.05 / n)))
q <- p.adjust(ok$p, "BH")
cat(sprintf("  min FDR = %.3f ; CpGs with FDR<0.10: %d\n", min(q), sum(q < 0.10)))

cat("\n=== power floor ===\n")
med_se <- median(ok$se)
zb <- qnorm(1 - 0.025 / n)
cat(sprintf("  median SE(beta) %.4f\n", med_se))
cat(sprintf("  80%% power at Bonferroni alpha needs |beta| >= %.3f (OR >= %.2f)\n",
            (zb + qnorm(0.8)) * med_se, exp((zb + qnorm(0.8)) * med_se)))
cat(sprintf("  observed median |beta| %.3f -> shortfall factor %.1fx\n",
            median(abs(ok$b)), (zb + qnorm(0.8)) * med_se / median(abs(ok$b))))

cat("\n=== locus-level (best CpG per gene) ===\n")
byg <- split(ok, ok$gene)
res <- do.call(rbind, lapply(names(byg), function(g) {
  x <- byg[[g]]
  data.frame(gene = g, nCpG = nrow(x), nSNP = length(unique(x$rsid)),
             p_min = min(x$p), cpg = x$cpg[which.min(x$p)],
             beta = x$b[which.min(x$p)], se = x$se[which.min(x$p)],
             F = x$F[which.min(x$p)])
}))
res$p_sidak <- 1 - (1 - res$p_min)^res$nSNP
res <- res[order(res$p_min), ]
print(format(res, digits = 3), row.names = FALSE)
write.table(res, "results/stage1b_locus_summary.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

writeLines(capture.output(sessionInfo()), "logs/sessioninfo_stage1b.txt")
