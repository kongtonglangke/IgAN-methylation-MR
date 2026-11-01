# DEFA1 colocalisation: does the result survive removal of the duplicated
# segment and of the lead variant's immediate neighbourhood?
#
# The concern being tested: cg00705004 lies in the 8p23.1 defensin cluster,
# where DEFA1/DEFA1B/DEFA3 form a 99.5%-identical 19 kb tandem repeat and the
# gene itself is copy-number variable.  If the mQTL were a copy-number readout
# rather than a local regulatory effect, PP.H4 should depend on variants inside
# the repeat, or collapse once the single lead SNP is removed.
#
# Three SNP sets are compared for each CpG:
#   all            - every shared SNP (baseline, reproduces the main analysis)
#   outside_dup    - drop SNPs inside the 99.5% repeat (hg19 6,839,960-6,878,169)
#   outside_dup_2k - as above, and additionally drop everything within 2 kb of
#                    the index SNP, so no single variant can carry the result

suppressMessages({
  library(coloc)
  library(jsonlite)
  library(data.table)
})

N_MQTL <- 27750
N_IGAN <- 822 + 499526
S_IGAN <- 822 / (822 + 499526)
P1 <- 1e-4; P2 <- 1e-4; P12 <- 1e-5
PALINDROME <- c("AT", "TA", "CG", "GC")

DUP_LO <- 6839960L     # hg19, genomicSuperDups fracMatch 0.9952
DUP_HI <- 6878169L
LEAD   <- 6822785L     # rs2738048, hg19

parse_pos <- function(snp) as.integer(sub("^chr[^:]+:([0-9]+):.*$", "\\1", snp))

align_one <- function(a1, a2, freq_m, ea, oa, eaf, beta_g) {
  a1 <- toupper(a1); a2 <- toupper(a2); ea <- toupper(ea); oa <- toupper(oa)
  if (any(is.na(c(a1, a2, ea, oa)))) return(NA_real_)
  if (paste0(a1, a2) %in% PALINDROME) {
    if (is.na(freq_m) || is.na(eaf)) return(NA_real_)
    if (abs(freq_m - eaf) < 0.10) return(beta_g)
    if (abs(freq_m - (1 - eaf)) < 0.10) return(-beta_g)
    return(NA_real_)
  }
  if (a1 == ea && a2 == oa) return(beta_g)
  if (a1 == oa && a2 == ea) return(-beta_g)
  NA_real_
}

build <- function(cpg, region) {
  mq <- as.data.table(fromJSON(file.path("data/coloc", sprintf("mqtl_%s.json", cpg))))
  ig <- fread(file.path("data/coloc", sprintf("finngen_%s.tsv", region)))
  mq <- mq[!is.na(rsid) & rsid != "" & !is.na(beta_a1) & !is.na(se) & se > 0]
  mq[, pos := parse_pos(snp)]
  mq <- mq[!duplicated(rsid)]
  mq2 <- mq[, .(rsid, a1, a2, beta_m = beta_a1, se_m = se,
                freq_m = freq_a1, p_m = pval, pos)]
  ig <- ig[rsid != "" & !is.na(beta) & !is.na(se) & se > 0]
  ig <- ig[!duplicated(rsid)]
  ig2 <- ig[, .(rsid, ea, oa, beta_g = beta, se_g = se, eaf,
                p_g = as.numeric(p), pos38)]
  m <- merge(mq2, ig2, by = "rsid")
  m[, beta_ig := mapply(align_one, a1, a2, freq_m, ea, oa, eaf, beta_g)]
  m <- m[!is.na(beta_ig)]
  m[, `:=`(eaf = suppressWarnings(as.numeric(eaf)),
           freq_m = suppressWarnings(as.numeric(freq_m)))]
  m <- m[is.finite(eaf) & is.finite(freq_m) & !is.na(pos) & !is.na(pos38)]
  m
}

run <- function(m, label) {
  d1 <- list(beta = as.numeric(m$beta_m), varbeta = as.numeric(m$se_m)^2,
             snp = as.character(m$rsid), position = as.integer(m$pos),
             type = "quant", N = N_MQTL, sdY = 1,
             MAF = pmin(m$freq_m, 1 - m$freq_m))
  d2 <- list(beta = as.numeric(m$beta_ig), varbeta = as.numeric(m$se_g)^2,
             snp = as.character(m$rsid), position = as.integer(m$pos),
             type = "cc", N = N_IGAN, s = S_IGAN,
             MAF = pmin(m$eaf, 1 - m$eaf))
  r <- tryCatch(suppressMessages(coloc.abf(d1, d2, p1 = P1, p2 = P2, p12 = P12)),
                error = function(e) NULL)
  if (is.null(r)) {
    cat(sprintf("  %-16s %6d  coloc failed\n", label, nrow(m)))
    return(invisible(NULL))
  }
  s <- r$summary
  pp <- as.data.table(r$results)
  setnames(pp, "SNP.PP.H4", "snp_pp", skip_absent=TRUE)
  pp <- pp[order(-snp_pp)]
  top <- pp[1]
  n_half <- sum(pp$snp_pp > 0.5)
  cat(sprintf("  %-16s %6d  H3=%.3f  H4=%.3f  | top=%s (%.2f)  SNPs w/ PP>0.5: %d\n",
              label, nrow(m), s[["PP.H3.abf"]], s[["PP.H4.abf"]],
              top$snp, top$snp_pp, n_half))
  invisible(list(h3 = s[["PP.H3.abf"]], h4 = s[["PP.H4.abf"]],
                 top = top$snp, top_pp = top$snp.PP.H4, n = nrow(m),
                 n_pp50 = n_half))
}

TARGETS <- list(
  c("cg00705004", "DEFA1_8p23",   "DEFA1 (H4 = 0.976)"),
  c("cg00539347", "HORMAD2_22q12", "HORMAD2 (H4 = 0.982)"),
  c("cg13682187", "CFH_CFHR_1q32", "CFH (H3-dominant)")
)

cat("Sensitivity of PP.H4 to removing the duplicated segment / the index SNP\n")
cat("duplicated block (hg19): ", DUP_LO, "-", DUP_HI, " fracMatch 0.9952\n\n")

rows <- list()
for (t in TARGETS) {
  cpg <- t[1]; region <- t[2]
  m0 <- build(cpg, region)
  cat(sprintf("%s  %s\n", cpg, t[3]))
  a <- run(m0, "all")
  m1 <- m0[!(pos >= DUP_LO & pos <= DUP_HI)]
  b <- run(m1, "outside_dup")
  m2 <- m1[abs(pos - LEAD) > 2000]
  cc <- run(m2, "outside_dup_2k")
  for (res in list(list("all", a), list("outside_dup", b), list("outside_dup_2k", cc))) {
    if (is.null(res[[2]])) next
    rows[[length(rows) + 1]] <- data.table(cpg = cpg, region = region,
                                           snpset = res[[1]], n = res[[2]]$n,
                                           H3 = res[[2]]$h3, H4 = res[[2]]$h4,
                                           top_snp = res[[2]]$top,
                                           top_pp = res[[2]]$top_pp,
                                           n_pp50 = res[[2]]$n_pp50)
  }
  cat("\n")
}

out <- rbindlist(rows)
fwrite(out, "results/coloc_cnv_sensitivity.tsv", sep = "\t")
cat("wrote results/coloc_cnv_sensitivity.tsv\n")
