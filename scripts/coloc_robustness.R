# Robustness of the two positive colocalisations.
#  1. sweep the p12 prior (how much belief in a shared variant is required)
#  2. report the SNP driving H4, i.e. the candidate shared causal variant
#  3. check whether that SNP is the IgAN lead variant for the region

suppressMessages({
  library(coloc); library(jsonlite); library(data.table)
})

N_MQTL <- 27750
N_IGAN <- 822 + 499526
S_IGAN <- 822 / (822 + 499526)
PAL <- c("AT", "TA", "CG", "GC")

CASES <- list(
  list(cpg = "cg00705004", region = "DEFA1_8p23",    lead = "rs2738048"),
  list(cpg = "cg00539347", region = "HORMAD2_22q12", lead = "rs12537"),
  list(cpg = "cg27665648", region = "HORMAD2_22q12", lead = "rs2412971")
)

align_one <- function(a1, a2, fm, ea, oa, eaf, bg) {
  a1 <- toupper(a1); a2 <- toupper(a2); ea <- toupper(ea); oa <- toupper(oa)
  if (any(is.na(c(a1, a2, ea, oa)))) return(NA_real_)
  if (paste0(a1, a2) %in% PAL) {
    if (is.na(fm) || is.na(eaf)) return(NA_real_)
    if (abs(fm - eaf) < 0.10) return(bg)
    if (abs(fm - (1 - eaf)) < 0.10) return(-bg)
    return(NA_real_)
  }
  if (a1 == ea && a2 == oa) return(bg)
  if (a1 == oa && a2 == ea) return(-bg)
  NA_real_
}

prep <- function(cpg, region) {
  mq <- as.data.table(fromJSON(sprintf("data/coloc/mqtl_%s.json", cpg)))
  ig <- fread(sprintf("data/coloc/finngen_%s.tsv", region))
  mq <- mq[!is.na(rsid) & rsid != "" & !is.na(beta_a1) & !is.na(se) & se > 0]
  mq[, pos := as.integer(sub("^chr[^:]+:([0-9]+):.*$", "\\1", snp))]
  mq <- mq[!duplicated(rsid)]
  ig <- ig[rsid != "" & !is.na(beta) & !is.na(se) & se > 0]
  ig <- ig[!duplicated(rsid)]
  m <- merge(mq[, .(rsid, a1, a2, beta_m = beta_a1, se_m = se,
                    freq_m = freq_a1, p_m = pval, pos)],
             ig[, .(rsid, ea, oa, beta_g = beta, se_g = se, eaf,
                    p_g = as.numeric(p))], by = "rsid")
  m[, beta_ig := mapply(align_one, a1, a2, freq_m, ea, oa, eaf, beta_g)]
  m <- m[!is.na(beta_ig)]
  m[, `:=`(eaf = suppressWarnings(as.numeric(eaf)),
           freq_m = suppressWarnings(as.numeric(freq_m)))]
  m <- m[is.finite(eaf) & is.finite(freq_m)]
  m
}

run <- function(m, p12) {
  d1 <- list(beta = as.numeric(m$beta_m), varbeta = as.numeric(m$se_m)^2,
             snp = as.character(m$rsid), position = as.integer(m$pos),
             type = "quant", N = N_MQTL, sdY = 1,
             MAF = pmin(m$freq_m, 1 - m$freq_m))
  d2 <- list(beta = as.numeric(m$beta_ig), varbeta = as.numeric(m$se_g)^2,
             snp = as.character(m$rsid), position = as.integer(m$pos),
             type = "cc", N = N_IGAN, s = S_IGAN,
             MAF = pmin(m$eaf, 1 - m$eaf))
  suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = p12))
}

cat(sprintf("%-12s %-16s %6s | %s\n", "cpg", "region", "nSNP",
            paste(sprintf("p12=%.0e", c(1e-6, 1e-5, 1e-4)), collapse = "  ")))
cat(strrep("-", 78), "\n")

for (cs in CASES) {
  m <- prep(cs$cpg, cs$region)
  if (nrow(m) < 20) { cat(sprintf("%-12s too few SNPs\n", cs$cpg)); next }
  h4s <- sapply(c(1e-6, 1e-5, 1e-4), function(p) run(m, p)$summary["PP.H4.abf"])
  cat(sprintf("%-12s %-16s %6d | %s\n", cs$cpg, cs$region, nrow(m),
              paste(sprintf("%.3f", h4s), collapse = "    ")))

  res <- run(m, 1e-5)
  r <- as.data.table(res$results)
  r <- r[order(-SNP.PP.H4)]
  top <- head(r, 5)
  cat("    top shared-variant candidates (SNP.PP.H4):\n")
  for (i in seq_len(nrow(top))) {
    s <- top$snp[i]
    is_lead <- s == cs$lead
    row <- m[rsid == s]
    cat(sprintf("      %-14s H4=%.3f  %s p_mqtl=%.1e p_igan=%.3f%s\n",
                s, top$SNP.PP.H4[i],
                ifelse(is_lead, "[IgAN lead]", "           "),
                row$p_m, row$p_g,
                ifelse(is_lead, "   <-- candidate causal SNP", "")))
  }
  lead_in <- m[rsid == cs$lead]
  if (nrow(lead_in)) {
    cat(sprintf("    IgAN lead %s: p_mqtl=%.1e p_igan=%.3g  rank_by_H4=%d / %d\n",
                cs$lead, lead_in$p_m[1], lead_in$p_g[1],
                which(r$snp == cs$lead), nrow(r)))
  } else {
    cat(sprintf("    IgAN lead %s NOT in the shared SNP set\n", cs$lead))
  }
  cat("\n")
}
