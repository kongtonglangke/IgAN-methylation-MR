# Prior sensitivity of PP.H4 for the two colocalized CpGs
# ---------------------------------------------------------------------------
# Re-uses the exact alignment and merge logic of coloc_igan_mqtl.R so that the
# default-prior cell reproduces the published PP.H4.  Anything else would make
# the sensitivity table incomparable with Table 2.

suppressMessages({
  library(coloc)
  library(jsonlite)
  library(data.table)
})

MQTL_DIR <- "data/coloc"
IGAN_DIR <- "data/coloc"
OUT      <- "results/coloc_prior_sensitivity.tsv"
N_MQTL   <- 27750
N_IGAN   <- 822 + 499526
S_IGAN   <- 822 / (822 + 499526)
PALINDROME <- c("AT", "TA", "CG", "GC")

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

prep <- function(cpg, region) {
  mq <- as.data.table(fromJSON(file.path(MQTL_DIR, sprintf("mqtl_%s.json", cpg))))
  ig <- fread(file.path(IGAN_DIR, sprintf("finngen_%s.tsv", region)))
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
           freq_m = suppressWarnings(as.numeric(freq_m)),
           pos = suppressWarnings(as.integer(pos)),
           pos38 = suppressWarnings(as.integer(pos38)))]
  m <- m[is.finite(eaf) & is.finite(freq_m) & !is.na(pos) & !is.na(pos38)]
  d1 <- list(beta = as.numeric(m$beta_m), varbeta = as.numeric(m$se_m)^2,
             snp = as.character(m$rsid), position = as.integer(m$pos),
             type = "quant", N = N_MQTL, sdY = 1,
             MAF = pmin(m$freq_m, 1 - m$freq_m))
  d2 <- list(beta = as.numeric(m$beta_ig), varbeta = as.numeric(m$se_g)^2,
             snp = as.character(m$rsid), position = as.integer(m$pos),
             type = "cc", N = N_IGAN, s = S_IGAN,
             MAF = pmin(m$eaf, 1 - m$eaf))
  list(d1 = d1, d2 = d2, n = nrow(m))
}

TARGETS <- list(c("cg00705004", "DEFA1_8p23"),
                c("cg00539347", "HORMAD2_22q12"))

P1V  <- c(1e-5, 1e-4, 1e-3)
P12V <- c(1e-6, 1e-5, 1e-4)

rows <- list()
prepped <- lapply(TARGETS, function(x) prep(x[1], x[2]))

for (i in seq_along(TARGETS)) {
  cpg <- TARGETS[[i]][1]; reg <- TARGETS[[i]][2]
  d <- prepped[[i]]
  for (p1 in P1V) for (p12 in P12V) {
    r <- suppressMessages(coloc.abf(d$d1, d$d2, p1 = p1, p2 = p1, p12 = p12))
    s <- r$summary
    rows[[length(rows) + 1]] <- data.table(
      cpg = cpg, region = reg, n_shared = d$n,
      p1p2 = p1, p12 = p12,
      PP.H0 = s[["PP.H0.abf"]], PP.H1 = s[["PP.H1.abf"]],
      PP.H2 = s[["PP.H2.abf"]], PP.H3 = s[["PP.H3.abf"]],
      PP.H4 = s[["PP.H4.abf"]])
  }
}
res <- rbindlist(rows)
fwrite(res, OUT, sep = "\t")

cat("== 默认先验格子（p1=p2=1e-4, p12=1e-5）==\n")
print(res[p1p2 == 1e-4 & p12 == 1e-5, .(cpg, n_shared, PP.H3, PP.H4)])
cat("\n== PP.H4 随先验的变化 ==\n")
print(dcast(res, cpg + p12 ~ p1p2, value.var = "PP.H4"))
cat("\n== 每个 CpG 的 PP.H4 范围 ==\n")
print(res[, .(min = min(PP.H4), max = max(PP.H4),
              below_0.8 = sum(PP.H4 < 0.8), n = .N), by = cpg])
