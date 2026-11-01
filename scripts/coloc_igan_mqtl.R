# Colocalisation of IgAN risk loci with target-gene methylation
# ---------------------------------------------------------------------------
# dataset 1: GoDMC cis-mQTL, beta = SD change in rank-normalised methylation
#            per allele -> type "quant", sdY = 1
# dataset 2: Sakaue 2021 IgAN GWAS (GCST90018866) -> type "cc"
#
# Allele alignment is done explicitly before coloc.abf: the mQTL beta is per
# copy of a1, the GWAS beta is per copy of effect_allele.  Palindromic SNPs are
# oriented by allele frequency and dropped when frequency is uninformative,
# because a sign error inflates or destroys PP.H4.

suppressMessages({
  library(coloc)
  library(jsonlite)
  library(data.table)
})

MQTL_DIR <- "data/coloc"
IGAN_DIR <- "data/coloc"
OUT      <- "results/coloc_summary.tsv"

N_MQTL <- 27750        # GoDMC meta-analysis sample size (README)

# GWAS side.  Two IgAN datasets were considered:
#   Sakaue 2021 GCST90018866  - the only IgAN GWAS with public summary
#       statistics, but its signal is MHC-dominated: none of the seven target
#       regions reaches p<5e-8 (minimum p 2e-4 to 5e-5), so coloc cannot work.
#   FinnGen R12 N14_IGA_NEPHROPATHY - registry/biopsy based, 822 cases and
#       499,526 controls, full summary statistics with a tabix index.  Used here.
GWAS_TAG <- "finngen"
N_IGAN <- 822 + 499526
S_IGAN <- 822 / (822 + 499526)

P1 <- 1e-4
P2 <- 1e-4
P12 <- 1e-5

PALINDROME <- c("AT", "TA", "CG", "GC")

parse_pos <- function(snp) as.integer(sub("^chr[^:]+:([0-9]+):.*$", "\\1", snp))

# Re-sign one IgAN beta onto the mQTL a1 allele.
# Returns NA when the alleles cannot be oriented.
align_one <- function(a1, a2, freq_m, ea, oa, eaf, beta_g) {
  a1 <- toupper(a1); a2 <- toupper(a2)
  ea <- toupper(ea); oa <- toupper(oa)
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

run_one <- function(cpg, region) {
  mf <- file.path(MQTL_DIR, sprintf("mqtl_%s.json", cpg))
  # one GWAS window per region: the CpGs of a region share the same query
  gf <- file.path(IGAN_DIR, sprintf("%s_%s.tsv", GWAS_TAG, region))
  if (!file.exists(mf) || !file.exists(gf)) {
    return(data.table(cpg = cpg, region = region, status = "missing_input"))
  }
  mq <- as.data.table(fromJSON(mf))
  ig <- fread(gf)
  if (nrow(mq) < 20 || nrow(ig) < 20) {
    return(data.table(cpg = cpg, region = region, status = "too_few_snps"))
  }

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
  if (nrow(m) < 20) {
    return(data.table(cpg = cpg, region = region,
                      status = "too_few_shared_snps", n_shared = nrow(m)))
  }

  m[, beta_ig := mapply(align_one, a1, a2, freq_m, ea, oa, eaf, beta_g)]
  n_drop <- sum(is.na(m$beta_ig))
  m <- m[!is.na(beta_ig)]
  # coloc refuses MAF/position vectors containing NA, and the GWAS columns are
  # read as character, so coerce and drop before building the datasets
  m[, `:=`(eaf = suppressWarnings(as.numeric(eaf)),
           freq_m = suppressWarnings(as.numeric(freq_m)),
           pos = suppressWarnings(as.integer(pos)),
           pos38 = suppressWarnings(as.integer(pos38)))]
  n_na <- sum(!is.finite(m$eaf) | !is.finite(m$freq_m) |
                is.na(m$pos) | is.na(m$pos38))
  m <- m[is.finite(eaf) & is.finite(freq_m) & !is.na(pos) & !is.na(pos38)]
  if (nrow(m) < 20) {
    return(data.table(cpg = cpg, region = region, status = "too_few_after_na",
                      n_shared = nrow(m) + n_na + n_drop, n_dropped = n_drop))
  }

  # Both datasets are aligned by rsID, so position is only auxiliary and must
  # be the *same* vector on both sides: d1 is GRCh37 (GoDMC) and d2 is GRCh38
  # (GWAS Catalog), and coloc rejects a position sequence that differs between
  # datasets even when the SNP names match.
  d1 <- list(beta = as.numeric(m$beta_m), varbeta = as.numeric(m$se_m)^2,
             snp = as.character(m$rsid), position = as.integer(m$pos),
             type = "quant", N = N_MQTL, sdY = 1,
             MAF = pmin(m$freq_m, 1 - m$freq_m))
  d2 <- list(beta = as.numeric(m$beta_ig), varbeta = as.numeric(m$se_g)^2,
             snp = as.character(m$rsid), position = as.integer(m$pos),
             type = "cc", N = N_IGAN, s = S_IGAN,
             MAF = pmin(m$eaf, 1 - m$eaf))

  res <- tryCatch(
    suppressMessages(coloc.abf(dataset1 = d1, dataset2 = d2,
                               p1 = P1, p2 = P2, p12 = P12)),
    error = function(e) {
      message(sprintf("[coloc error] %s: %s", cpg, conditionMessage(e)))
      NULL
    })
  if (is.null(res)) {
    return(data.table(cpg = cpg, region = region, status = "coloc_error",
                      n_shared = nrow(m)))
  }
  s <- as.list(res$summary)
  sens <- tryCatch({
    r <- sensitivity(res, "H4 > 0.8")
    paste(sprintf("%s:%.2f", rownames(r$sensitivity), r$sensitivity[, 1]),
          collapse = "; ")
  }, error = function(e) NA_character_)

  data.table(cpg = cpg, region = region, status = "ok",
             n_shared = nrow(m), n_dropped = n_drop,
             n_mqtl_sig = sum(m$p_m < 5e-8),
             n_igan_sig = sum(m$p_g < 5e-8),
             PP.H0 = s$PP.H0.abf, PP.H1 = s$PP.H1.abf, PP.H2 = s$PP.H2.abf,
             PP.H3 = s$PP.H3.abf, PP.H4 = s$PP.H4.abf,
             nsnps = s$nsnps, sensitivity = sens)
}

TARGETS <- list(
  c("cg13682187", "CFH_CFHR_1q32"), c("cg11348583", "CFH_CFHR_1q32"),
  c("cg13441112", "CFH_CFHR_1q32"), c("cg07209298", "CFH_CFHR_1q32"),
  c("cg00705004", "DEFA1_8p23"),    c("cg07297068", "DEFA1_8p23"),
  c("cg19304352", "DEFA1_8p23"),    c("cg05969298", "DEFA1_8p23"),
  c("cg03100639", "ITGAM_ITGAX_16p11"), c("cg10451425", "ITGAM_ITGAX_16p11"),
  c("cg02846316", "ITGAM_ITGAX_16p11"), c("cg04018474", "ITGAM_ITGAX_16p11"),
  c("cg27665648", "HORMAD2_22q12"), c("cg00539347", "HORMAD2_22q12"),
  c("cg23386659", "MADCAM1_19p13"), c("cg27535410", "MADCAM1_19p13"),
  c("cg17448192", "ST6GAL1_3q27"),  c("cg25599673", "ST6GAL1_3q27"),
  c("cg02671941", "TNFRSF13B_17p11")
)

out <- rbindlist(lapply(TARGETS, function(t) run_one(t[1], t[2])), fill = TRUE)
if (!"PP.H4" %in% names(out)) out[, PP.H4 := NA_real_]
setorder(out, region, -PP.H4, na.last = TRUE)
fwrite(out, OUT, sep = "\t")
cat(sprintf("\nwrote %s  (%d CpGs)\n\n", OUT, nrow(out)))

ok <- out[status == "ok"]
if (!nrow(ok)) {
  cat("no successful colocalisations\n")
  print(out[, .(cpg, region, status)])
} else {
  cat(sprintf("%-12s %-20s %6s %6s %6s %6s %6s %6s  %s\n",
              "cpg", "region", "H0", "H1", "H2", "H3", "H4", "nSNP", "tol"))
  for (i in seq_len(nrow(ok))) {
    r <- ok[i]
    cat(sprintf("%-12s %-20s %6.3f %6.3f %6.3f %6.3f %6.3f %6d  %s\n",
                r$cpg, r$region, r$PP.H0, r$PP.H1, r$PP.H2, r$PP.H3,
                r$PP.H4, r$nsnps, ifelse(is.na(r$n_dropped), "-", r$n_dropped)))
  }
  cat(sprintf("\nPP.H4 > 0.8 : %d / %d\n", sum(ok$PP.H4 > 0.8), nrow(ok)))
  cat(sprintf("PP.H4 > 0.5 : %d / %d\n", sum(ok$PP.H4 > 0.5), nrow(ok)))
  cat(sprintf("PP.H3 > 0.8 : %d / %d  (distinct causal variants)\n",
              sum(ok$PP.H3 > 0.8), nrow(ok)))
}

bad <- out[status != "ok"]
if (nrow(bad)) {
  cat("\n--- not analysed ---\n")
  for (i in seq_len(nrow(bad))) {
    cat(sprintf("  %-12s %-20s %-22s shared=%s\n", bad$cpg[i], bad$region[i],
                bad$status[i],
                if ("n_shared" %in% names(bad)) bad$n_shared[i] else "-"))
  }
}
