# 09_sensitivity.R
# PRE-REGISTERED (written 6 Oct 2026, before running):
#  Placebo is NON-NULL if any ring's 95% CI excludes 0 with |ATT| >= 0.25 C, or if
#    >= 2 rings' CIs exclude 0 in the same direction. Non-null => stop the analysis.
#  Negative control FAILS if the 95% CI of the 1975 built-up difference excludes 0
#    and |standardized difference| >= 0.1.
# Parts: A negative control, B Rosenbaum bounds, C far-ring diagnostic (extra, not in
# original spec), D placebo parks. All on the headline run (built neighbourhoods,
# urban parks, exogenous covariates). Runnable on its own.
source("code/_common.R")   # run from the repository folder
suppressPackageStartupMessages({ library(MatchIt); library(cobalt); library(fixest); library(future.apply); library(rbounds) })
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1"); setFixest_nthreads(1)
plan(multisession, workers = 6)
set.seed(SEED)
pool <- readRDS(file.path(DIR_PROC, "pool_locations.rds")); long <- readRDS(file.path(DIR_PROC, "pool_lst_long.rds"))
RINGS <- c(0, 50, 100, 150, 200, 300, 400, 500)
RING_LAB <- c(paste0(RINGS[-8], "-", RINGS[-1], " m"), "park interior")
COV_A <- c("elev_m", "slope_deg", "northness", "eastness", "dist_water_m", "dist_cbd_m", "dist_coast_m")
BLDG_MIN <- 0.02; NBOOT <- 500L
parks <- st_read(file.path(DIR_PROC, "parks_main_urban.gpkg"), quiet = TRUE)
yr <- data.table(park_id = as.character(parks$UNIT_ID), yr_est = parks$yr_est)
b75 <- setNames(pool$built1975, pool$cell)
out <- list()
bootblocks <- function(v, blk, B = NBOOT) { bl <- unique(blk); idx <- split(seq_along(v), blk)
  vapply(seq_len(B), function(b) { set.seed(SEED + b); mean(v[unlist(idx[sample(bl, length(bl), TRUE)], use.names = FALSE)]) }, 0) }
pairs_of <- function(md) { md <- copy(md); md[, park_id := as.character(park_id)]
  merge(md[treat == 1, .(subclass, scene, cell_t = cell, park_id, block, lst_t = lst)],
        md[treat == 0, .(subclass, cell_c = cell, lst_c = lst)], by = "subclass") }

## A. Negative control: 1975 built-up, pairs whose park was established after 1975
sd75 <- sd(pool[!is.na(ring_main) & bldg_dens150 >= BLDG_MIN, built1975])
for (k in 1:8) {
  pr <- merge(pairs_of(readRDS(file.path(DIR_PROC, "matched", paste0("built_main_exog_ring", k, ".rds")))), yr, by = "park_id")
  pr <- pr[!is.na(yr_est) & yr_est > 1975]
  pr[, d75 := b75[as.character(cell_t)] - b75[as.character(cell_c)]]
  loc <- pr[, .(d = mean(d75)), by = .(cell_t, block)]
  if (nrow(loc) < 20) { logmsg("NegCtrl ", RING_LAB[k], ": too few post-1975 pairs (", nrow(loc), ")"); next }
  bs <- bootblocks(loc$d, loc$block); ci <- quantile(bs, c(.025, .975))
  est <- mean(loc$d); smd <- est / sd75
  fail <- (ci[1] > 0 | ci[2] < 0) & abs(smd) >= 0.1
  out$neg[[k]] <- data.table(ring = RING_LAB[k], treated_locations = nrow(loc), parks = uniqueN(pr$park_id),
                             diff_built1975_pp = 100 * est, ci_lo_pp = 100 * ci[1], ci_hi_pp = 100 * ci[2], std_diff = smd, verdict = ifelse(fail, "FAIL", "pass"))
  logmsg(sprintf("NegCtrl %-14s | %d locs, %d parks | diff %.2f pp [%.2f, %.2f] | std diff %.3f | %s",
                 RING_LAB[k], nrow(loc), uniqueN(pr$park_id), 100 * est, 100 * ci[1], 100 * ci[2], smd, ifelse(fail, "FAIL", "pass")))
}
fwrite(rbindlist(out$neg), file.path(DIR_OUT, "sensitivity_negative_control.csv"))

## B. Rosenbaum bounds (Wilcoxon signed rank), pairs aggregated to one value per treated location
for (k in 1:8) {
  pr <- pairs_of(readRDS(file.path(DIR_PROC, "matched", paste0("built_main_exog_ring", k, ".rds"))))
  loc <- pr[, .(t = mean(lst_t), c = mean(lst_c)), by = cell_t]
  ps <- psens(x = loc$c, y = loc$t, Gamma = 6, GammaInc = 0.05)$bounds
  g  <- ps$Gamma[which(ps[["Upper bound"]] > 0.05)[1]]
  out$ros[[k]] <- data.table(ring = RING_LAB[k], treated_locations = nrow(loc), gamma_critical = ifelse(is.na(g), ">6", as.character(g)))
  logmsg(sprintf("Rosenbaum %-14s | %d locs | result stays significant (p<0.05) up to Gamma = %s", RING_LAB[k], nrow(loc), ifelse(is.na(g), ">6", sprintf("%.2f", g - 0.05))))
}
fwrite(rbindlist(out$ros), file.path(DIR_OUT, "sensitivity_rosenbaum.csv"))

## shared estimator (same rules as stage 07)
att_diff <- function(d) d[, weighted.mean(lst[treat == 1], weights[treat == 1]) - weighted.mean(lst[treat == 0], weights[treat == 0])]
att_dr <- function(d, covs) coef(feols(as.formula(paste("lst ~ treat +", paste(covs, collapse = " + "), "| scene")), data = d, weights = ~weights, notes = FALSE, warn = FALSE))[["treat"]]
boot <- function(d, covs) { bl <- unique(d$block); idx <- split(seq_len(nrow(d)), d$block)
  future_sapply(SEED + seq_len(NBOOT), function(s) { fixest::setFixest_nthreads(1); set.seed(s); db <- d[unlist(idx[sample(bl, length(bl), TRUE)], use.names = FALSE)]
    if (length(unique(db$treat)) < 2) return(c(NA, NA)); c(att_diff(db), tryCatch(att_dr(db, covs), error = function(e) NA)) },
    future.seed = TRUE, future.packages = c("data.table", "fixest")) }
est_ring <- function(Tk, Ck, covs, tag, lab) {
  loc <- rbind(Tk[, c("cell", "block", covs), with = FALSE][, treat := 1L], Ck[, c("cell", "block", covs), with = FALSE][, treat := 0L])
  loc[, ps_logit := predict(glm(as.formula(paste("treat ~", paste(covs, collapse = " + "))), family = binomial, data = loc), type = "link")]
  cr <- range(loc[treat == 0, ps_logit]); ntr <- loc[treat == 1 & (ps_logit < cr[1] | ps_logit > cr[2]), .N]
  loc <- loc[!(treat == 1 & (ps_logit < cr[1] | ps_logit > cr[2]))]
  d <- merge(long, loc, by = "cell"); d[, scene := factor(scene)]; nT <- sum(d$treat == 1)
  f <- as.formula(paste("treat ~", paste(covs, collapse = " + ")))
  m <- matchit(update(f, . ~ . + scene), data = d, method = "cem", estimand = "ATT"); meth <- "CEM"
  kept <- sum(m$weights[d$treat == 1] > 0) / nT
  if (kept < 0.40) { m <- matchit(f, data = d, method = "nearest", distance = "glm", link = "linear.logit", mahvars = as.formula(paste("~", paste(covs, collapse = " + "))),
                                  caliper = 0.2, std.caliper = TRUE, exact = ~scene, replace = FALSE, ratio = 1); meth <- "NN"; kept <- sum(m$weights[d$treat == 1] > 0) / nT }
  bt <- bal.tab(m, stats = c("m", "v"), un = TRUE, s.d.denom = "treated")$Balance; bt <- bt[rownames(bt) %in% covs, ]
  md <- as.data.table(match.data(m)); a1 <- att_diff(md); a2 <- att_dr(md, covs)
  bs <- boot(md[, c("lst", "treat", "weights", "scene", "block", covs), with = FALSE], covs)
  ci <- quantile(bs[1, ], c(.025, .975), na.rm = TRUE)
  unb <- rownames(bt)[abs(bt$Diff.Adj) >= .1 | bt$V.Ratio.Adj < .5 | bt$V.Ratio.Adj > 2]
  r <- data.table(analysis = tag, ring = lab, method = meth, treated_locations = nrow(Tk), trimmed = ntr, share_retained = round(kept, 3),
                  att = a1, se = sd(bs[1, ], na.rm = TRUE), ci_lo = ci[1], ci_hi = ci[2], att_dr = a2, se_dr = sd(bs[2, ], na.rm = TRUE),
                  max_abs_smd = round(max(abs(bt$Diff.Adj)), 3), unbalanced = ifelse(length(unb), paste(unb, collapse = "; "), "none"))
  logmsg(sprintf("%s %-14s | %s | n_T %d | kept %.0f%% | ATT %.2f C [%.2f, %.2f] | DR %.2f | max|SMD| %.3f | unbalanced: %s",
                 tag, lab, meth, nrow(Tk), 100 * kept, a1, ci[1], ci[2], a2, max(abs(bt$Diff.Adj)), r$unbalanced)); r
}

## C. Far-ring diagnostic (extra): rings 0-400 m vs the 400-500 m ring as control
B <- pool[!is.na(ring_main) & bldg_dens150 >= BLDG_MIN]
for (k in 1:6) out$far[[k]] <- est_ring(B[ring_main == k], B[ring_main == 7], COV_A, "FarRing", RING_LAB[k])
fwrite(rbindlist(out$far), file.path(DIR_OUT, "diagnostic_vs_400_500m_ring.csv"))

## D. Placebo parks: real urban-park shapes, random position + rotation, wholly inside control territory
set.seed(SEED + 1)
tm <- rast(file.path(DIR_PROC, "grid_template.tif")); uav <- vect(st_read(file.path(DIR_PROC, "study_area.gpkg"), quiet = TRUE))
zone <- (rast(file.path(DIR_PROC, "cov", "dist_any_green.tif")) > 800) & (rast(file.path(DIR_PROC, "cov", "dist_water_m.tif")) > 800) &
        (rast(file.path(DIR_PROC, "cov", "bldg_dens150.tif")) >= BLDG_MIN) & (rast(file.path(DIR_PROC, "cov", "perm_water.tif")) == 0)
zone <- mask(classify(zone, cbind(0, NA)), uav); zc <- cells(zone)
cp <- pool[control == TRUE & bldg_dens150 >= BLDG_MIN]
pts <- st_as_sf(cp[, .(x, y)], coords = c("x", "y"), crs = CRS_WORK)
shapes <- st_geometry(parks); covered <- rep(FALSE, nrow(pts)); placed <- list(); src <- c(); TARGET <- 0.40
for (i in sample(length(shapes))) {
  if (mean(covered) >= TARGET) break
  g <- shapes[i]; cxy <- as.numeric(st_coordinates(st_centroid(g)))
  for (a in 1:30) {
    p <- xyFromCell(zone, zc[sample.int(length(zc), 1)]); th <- runif(1, 0, 2 * pi)
    g2 <- st_set_crs((st_set_crs(g, NA) - cxy) * matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2, 2) + as.numeric(p), CRS_WORK)
    v <- terra::extract(zone, vect(g2), touches = TRUE, ID = FALSE)[[1]]
    if (!length(v) || anyNA(v)) next
    if (length(placed) && any(lengths(st_intersects(g2, do.call(c, placed))) > 0)) next
    placed[[length(placed) + 1]] <- g2; src <- c(src, i)
    covered <- covered | as.vector(st_is_within_distance(pts, g2, 500, sparse = FALSE)); break
  }
}
pl <- st_sf(src_unit = parks$UNIT_ID[src], area_ha = as.numeric(st_area(do.call(c, placed))) / 1e4, geometry = do.call(c, placed))
st_write(pl, file.path(DIR_PROC, "placebo_parks.gpkg"), delete_dsn = TRUE, quiet = TRUE)
logmsg("Placebo parks placed: ", nrow(pl), " | control locations within 500 m: ", round(100 * mean(covered), 1), "%")
cat("Placebo vs real urban-park area (ha) quantiles:\n"); print(rbind(placebo = round(quantile(pl$area_ha, c(.1, .25, .5, .75, .9)), 2), real = round(quantile(parks$area_ha, c(.1, .25, .5, .75, .9)), 2)))
cp[, d_pl := as.numeric(st_distance(pts, st_union(pl)))]
cp[, ring_pl := as.integer(cut(d_pl, RINGS, right = TRUE))]
Cpl <- cp[d_pl > 800]; logmsg("Placebo controls (>800 m from any placebo park): ", nrow(Cpl))
for (k in 1:8) {
  Tk <- if (k <= 7) cp[!is.na(ring_pl) & d_pl > 0 & ring_pl == k] else cp[d_pl == 0]
  if (nrow(Tk) < 30) { logmsg("Placebo ", RING_LAB[k], ": only ", nrow(Tk), " treated locations, skipped"); next }
  out$pla[[k]] <- est_ring(Tk, Cpl, COV_A, "Placebo", RING_LAB[k])
}
P <- rbindlist(out$pla); P[, excl0 := ci_lo > 0 | ci_hi < 0]; P[, ci_halfwidth := (ci_hi - ci_lo) / 2]
fwrite(P, file.path(DIR_OUT, "sensitivity_placebo.csv"))
nonnull <- any(P$excl0 & abs(P$att) >= 0.25) | any(table(sign(P$att[P$excl0])) >= 2)
logmsg("Placebo median 95% CI half-width: ", round(median(P$ci_halfwidth), 2), " C (smallest effect the test could reliably flag is roughly this size)")
logmsg("PLACEBO VERDICT: ", ifelse(nonnull, "NON-NULL -> STOP per pre-registered rule", "null (passes)"))
logmsg("STAGE_09_DONE")
