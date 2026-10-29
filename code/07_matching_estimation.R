# 07_matching_estimation.R
# Per run (park definition x covariate set) and per ring: common-support trim,
# CEM (exact on scene) with Mahalanobis-NN fallback, balance diagnostics,
# matched-difference and doubly robust ATT, 1 km block bootstrap SEs.
# Runnable on its own: Rscript 07_matching_estimation.R [run_label]
source("code/_common.R")   # run from the repository folder
suppressPackageStartupMessages({ library(MatchIt); library(cobalt); library(fixest); library(future.apply) })
set.seed(SEED)
args <- commandArgs(trailingOnly = TRUE)
pool <- readRDS(file.path(DIR_PROC, "pool_locations.rds"))
long <- readRDS(file.path(DIR_PROC, "pool_lst_long.rds"))
RINGS <- c(0, 50, 100, 150, 200, 300, 400, 500)
RING_LAB <- c(paste0(RINGS[-8], "-", RINGS[-1], " m"), "park interior")
COV_A <- c("elev_m", "slope_deg", "northness", "eastness", "dist_water_m", "dist_cbd_m", "dist_coast_m")
COV_B <- c(COV_A, "bldg_dens150")
BLDG_MIN <- 0.02   # pre-registered 6 Oct 2026, before any ATT was estimated: built-neighbourhood scope
RUNS <- list(
  built_main_exog        = list(ring = "ring_main", int = "interior_main", park = "park_main", covs = COV_A, restrict = TRUE),
  built_main_plus_bldg   = list(ring = "ring_main", int = "interior_main", park = "park_main", covs = COV_B, restrict = TRUE),
  built_allgreen_exog    = list(ring = "ring_all",  int = "interior_all",  park = "park_all",  covs = COV_A, restrict = TRUE),
  unrestr_main_exog      = list(ring = "ring_main", int = "interior_main", park = "park_main", covs = COV_A, restrict = FALSE),
  unrestr_main_plus_bldg = list(ring = "ring_main", int = "interior_main", park = "park_main", covs = COV_B, restrict = FALSE))
if (length(args)) RUNS <- RUNS[args]
TARGET_T <- 40000L; TARGET_C <- 150000L; NBOOT <- 500L
dir.create(file.path(DIR_PROC, "matched"), showWarnings = FALSE)
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")   # inherited by bootstrap workers
setFixest_nthreads(1)
plan(multisession, workers = 6)   # 6 single-threaded workers: fits in laptop RAM

## pre-registered subsampling (locations, after thinning) if pools exceed targets
ctrl_all <- pool[control == TRUE]
if (nrow(ctrl_all) > TARGET_C) { ctrl_all <- ctrl_all[sample(.N, TARGET_C)]; logmsg("Controls subsampled to ", TARGET_C) }

att_diff <- function(d) d[, weighted.mean(lst[treat == 1], weights[treat == 1]) - weighted.mean(lst[treat == 0], weights[treat == 0])]
att_dr   <- function(d, covs) coef(feols(as.formula(paste("lst ~ treat +", paste(covs, collapse = " + "), "| scene")),
                                         data = d, weights = ~weights, notes = FALSE, warn = FALSE))[["treat"]]
boot <- function(d, covs, B = NBOOT) {
  bl <- unique(d$block); idx <- split(seq_len(nrow(d)), d$block)
  seeds <- SEED + seq_len(B)
  r <- future_sapply(seeds, function(s) {
    fixest::setFixest_nthreads(1); set.seed(s); pick <- sample(bl, length(bl), replace = TRUE)
    db <- d[unlist(idx[pick], use.names = FALSE)]
    if (length(unique(db$treat)) < 2) return(c(NA, NA))
    c(att_diff(db), tryCatch(att_dr(db, covs), error = function(e) NA))
  }, future.seed = TRUE, future.packages = c("data.table", "fixest"))
  r
}

res <- list(); bal <- list(); trim <- list()
for (rn in names(RUNS)) {
  R <- RUNS[[rn]]; covs <- R$covs
  trt_all  <- pool[(!is.na(get(R$ring)) & (!R$restrict | bldg_dens150 >= BLDG_MIN)) | get(R$int) == TRUE]
  ctrl_run <- if (R$restrict) ctrl_all[bldg_dens150 >= BLDG_MIN] else ctrl_all
  logmsg(rn, ": treated locations ", nrow(trt_all), " | control locations ", nrow(ctrl_run))
  if (nrow(trt_all) > TARGET_T) { trt_all <- trt_all[sample(.N, TARGET_T)]; logmsg(rn, ": treated locations subsampled to ", TARGET_T) }
  for (k in 1:8) {
    lab <- RING_LAB[k]
    Tk <- if (k <= 7) trt_all[get(R$ring) == k] else trt_all[get(R$int) == TRUE]
    if (nrow(Tk) < 30) { logmsg(rn, " | ", lab, ": only ", nrow(Tk), " treated locations, skipped"); next }
    loc <- rbind(Tk[, c("cell", "block", R$park, covs), with = FALSE][, treat := 1L],
                 ctrl_run[, c("cell", "block", R$park, covs), with = FALSE][, treat := 0L])
    setnames(loc, R$park, "park_id")
    ## common support on the propensity logit (location level)
    fps <- glm(as.formula(paste("treat ~", paste(covs, collapse = " + "))), family = binomial, data = loc)
    loc[, ps_logit := predict(fps, type = "link")]
    cr <- range(loc[treat == 0, ps_logit])
    loc[, trim_reason := fifelse(treat == 1 & ps_logit < cr[1], "propensity below control range",
                          fifelse(treat == 1 & ps_logit > cr[2], "propensity above control range", NA_character_))]
    trim[[length(trim) + 1]] <- loc[!is.na(trim_reason), .(run = rn, ring = lab, cell, trim_reason)]
    n_loc_T <- sum(loc$treat == 1); n_trim <- sum(!is.na(loc$trim_reason))
    loc <- loc[is.na(trim_reason)]
    d <- merge(long, loc, by = "cell"); d[, scene := factor(scene)]
    n_obs_T <- sum(d$treat == 1)
    f <- as.formula(paste("treat ~", paste(covs, collapse = " + ")))
    m <- matchit(update(f, . ~ . + scene), data = d, method = "cem", estimand = "ATT")   # factor scene => matched exactly by CEM
    method <- "CEM"
    kept_T <- sum(m$weights[d$treat == 1] > 0) / n_obs_T
    if (kept_T < 0.40) {
      logmsg(rn, " | ", lab, ": CEM kept ", round(100 * kept_T, 1), "% of treated -> Mahalanobis NN fallback")
      m <- matchit(f, data = d, method = "nearest", distance = "glm", link = "linear.logit",
                   mahvars = as.formula(paste("~", paste(covs, collapse = " + "))), caliper = 0.2, std.caliper = TRUE, exact = ~scene, replace = FALSE, ratio = 1)
      method <- "Mahalanobis NN (0.2 SD logit caliper, 1:1, no replacement)"
      kept_T <- sum(m$weights[d$treat == 1] > 0) / n_obs_T
    }
    bt <- bal.tab(m, stats = c("m", "v"), un = TRUE, binary = "std", s.d.denom = "treated")$Balance
    bt <- data.table(covariate = rownames(bt), bt)[covariate %in% covs]
    bt <- bt[, .(run = rn, ring = lab, covariate, smd_before = Diff.Un, smd_after = Diff.Adj, vr_before = V.Ratio.Un, vr_after = V.Ratio.Adj)]
    bal[[length(bal) + 1]] <- bt
    md <- as.data.table(match.data(m))
    saveRDS(md, file.path(DIR_PROC, "matched", paste0(rn, "_ring", k, ".rds")))
    a1 <- att_diff(md); a2 <- att_dr(md, covs)
    bs <- boot(md[, c("lst", "treat", "weights", "scene", "block", covs), with = FALSE], covs)
    se1 <- sd(bs[1, ], na.rm = TRUE); se2 <- sd(bs[2, ], na.rm = TRUE)
    ci1 <- quantile(bs[1, ], c(.025, .975), na.rm = TRUE); ci2 <- quantile(bs[2, ], c(.025, .975), na.rm = TRUE)
    unbal <- bt[abs(smd_after) >= 0.1 | vr_after < 0.5 | vr_after > 2, paste(covariate, collapse = "; ")]
    res[[length(res) + 1]] <- data.table(run = rn, ring = lab, method,
      treated_locations = n_loc_T, trimmed_locations = n_trim, treated_obs_after_trim = n_obs_T,
      treated_obs_matched = sum(md$treat == 1), share_treated_retained = round(kept_T, 3),
      control_obs_used = sum(md$treat == 0), blocks = uniqueN(md$block),
      att_matched_diff = a1, se_matched_diff = se1, ci95_lo_diff = ci1[1], ci95_hi_diff = ci1[2],
      att_doubly_robust = a2, se_dr = se2, ci95_lo_dr = ci2[1], ci95_hi_dr = ci2[2],
      boot_reps_ok = sum(!is.na(bs[1, ])), max_abs_smd_after = round(max(abs(bt$smd_after)), 3),
      unbalanced_covariates = ifelse(unbal == "", "none", unbal))
    logmsg(sprintf("%s | %-14s | %s | kept %.0f%% | ATT diff %.2f C (SE %.2f) | DR %.2f C (SE %.2f) | max|SMD| %.3f | unbalanced: %s",
                   rn, lab, ifelse(method == "CEM", "CEM", "NN"), 100 * kept_T, a1, se1, a2, se2, max(abs(bt$smd_after)), ifelse(unbal == "", "none", unbal)))
    fwrite(rbindlist(res), file.path(DIR_OUT, "results_att_by_ring_partial.csv")); fwrite(rbindlist(bal), file.path(DIR_OUT, "balance_table_partial.csv"))
    fwrite(rbindlist(trim), file.path(DIR_PROC, "trimmed_pixels_partial.csv"))
  }
}
R <- rbindlist(res); B <- rbindlist(bal); T <- rbindlist(trim)
fwrite(R, file.path(DIR_OUT, "results_att_by_ring.csv")); fwrite(B, file.path(DIR_OUT, "balance_table.csv"))
fwrite(T, file.path(DIR_PROC, "trimmed_pixels.csv"))
logmsg("Seed: ", SEED, " | bootstrap reps: ", NBOOT)
logmsg("STAGE_07_DONE")
