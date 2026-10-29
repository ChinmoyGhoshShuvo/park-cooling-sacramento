# 08_heterogeneity.R
# PRE-REGISTERED (6 Oct 2026, before running). Headline design: each ring vs matched
# pixels in the 400-500 m ring (built neighbourhoods, urban parks, exogenous covariates).
#  Bands: 0-50 m, 50-150 m ("beyond 50 m"), 150-200 m.
#  Park cooling distance = outer edge of the last contiguous band (from the edge) whose 95% CI < 0.
#  Area threshold = smallest area-class lower bound above which every class shows 50-150 m
#    cooling with 95% CI < 0. Classes (ha): 0.2-0.5, 0.5-1, 1-2, 2-4, 4-8, 8-16, >16.
#  Meta-regression: per-park 50-150 m ATT ~ log(area) + P/A + canopy share + internal water,
#    parks with >= 5 locations in the band, weights = locations, heteroskedasticity-robust SEs.
#  Canopy = ESA WorldCover 2021 tree class (10); water = class 80 (>= 0.05 ha counts as water).
source("code/_common.R")   # run from the repository folder
suppressPackageStartupMessages({ library(MatchIt); library(cobalt); library(fixest); library(future.apply); library(rstac); library(exactextractr) })
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1"); setFixest_nthreads(1); plan(multisession, workers = 6)
set.seed(SEED)
pool <- readRDS(file.path(DIR_PROC, "pool_locations.rds")); long <- readRDS(file.path(DIR_PROC, "pool_lst_long.rds"))
COV_A <- c("elev_m", "slope_deg", "northness", "eastness", "dist_water_m", "dist_cbd_m", "dist_coast_m")
BLDG_MIN <- 0.02; NBOOT <- 500L
parks <- st_read(file.path(DIR_PROC, "parks_main_urban.gpkg"), quiet = TRUE)

## 1. Park attributes
parks$perim_m  <- as.numeric(st_length(st_boundary(parks)))
parks$pa_ratio <- parks$perim_m / (parks$area_ha * 1e4)
bb <- as.numeric(st_bbox(st_transform(parks, 4326)))
s  <- stac("https://planetarycomputer.microsoft.com/api/stac/v1")
wit <- s |> stac_search(collections = "esa-worldcover", bbox = bb, datetime = "2021-01-01T00:00:00Z/2021-12-31T23:59:59Z") |> post_request() |> items_sign(sign_fn = sign_planetary_computer())
logmsg("WorldCover 2021 tiles: ", length(wit$features))
wc <- vrt(paste0("/vsicurl/", vapply(wit$features, function(f) f$assets$map$href, "")))
wc <- crop(wc, ext(bb[c(1, 3, 2, 4)]) + 0.01)
wc <- project(wc, CRS_WORK, method = "near", res = 10)
pv <- st_transform(parks, crs(wc))
parks$canopy_frac <- exact_extract(wc == 10, pv, "mean", progress = FALSE)
parks$water_ha    <- exact_extract(wc == 80, pv, "sum", progress = FALSE) * 100 / 1e4
parks$has_water   <- parks$water_ha >= 0.05
logmsg("Canopy share in parks: p10/p50/p90 = ", paste(round(quantile(parks$canopy_frac, c(.1, .5, .9), na.rm = TRUE), 3), collapse = "/"),
       " | parks with internal water: ", sum(parks$has_water, na.rm = TRUE), " of ", nrow(parks))
stopif(any(is.na(parks$canopy_frac)), "canopy share missing for some parks")
stopif(median(parks$canopy_frac) <= 0 | median(parks$canopy_frac) >= 0.9, "implausible park canopy share")

## 2. Far-ring matching (rings 1-4 and park interior vs the 400-500 m ring), keeping pairs
B  <- pool[!is.na(ring_main) & bldg_dens150 >= BLDG_MIN]; C7 <- B[ring_main == 7]
match_md <- function(Tk, Ck, covs) {
  loc <- rbind(Tk[, c("cell", "block", "park_main", covs), with = FALSE][, treat := 1L], Ck[, c("cell", "block", "park_main", covs), with = FALSE][, treat := 0L])
  loc[, ps_logit := predict(glm(as.formula(paste("treat ~", paste(covs, collapse = " + "))), family = binomial, data = loc), type = "link")]
  cr <- range(loc[treat == 0, ps_logit]); ntr <- loc[treat == 1 & (ps_logit < cr[1] | ps_logit > cr[2]), .N]
  loc <- loc[!(treat == 1 & (ps_logit < cr[1] | ps_logit > cr[2]))]
  d <- merge(long, loc, by = "cell"); d[, scene := factor(scene)]
  f <- as.formula(paste("treat ~", paste(covs, collapse = " + ")))
  m <- matchit(update(f, . ~ . + scene), data = d, method = "cem", estimand = "ATT"); meth <- "CEM"
  if (sum(m$weights[d$treat == 1] > 0) / sum(d$treat == 1) < 0.40) {
    m <- matchit(f, data = d, method = "nearest", distance = "glm", link = "linear.logit", mahvars = as.formula(paste("~", paste(covs, collapse = " + "))),
                 caliper = 0.2, std.caliper = TRUE, exact = ~scene, replace = FALSE, ratio = 1); meth <- "NN" }
  list(md = as.data.table(match.data(m)), m = m, trimmed = ntr, method = meth)
}
pairs_of <- function(md) merge(md[treat == 1, .(subclass, cell_t = cell, park_id = as.character(park_main), block, lst_t = lst)],
                               md[treat == 0, .(subclass, lst_c = lst)], by = "subclass")
ring_lab <- c("0-50 m", "50-100 m", "100-150 m", "150-200 m")
P <- list(); bal <- list()
for (k in 1:5) {
  Tk <- if (k <= 4) B[ring_main == k] else pool[interior_main == TRUE]
  lab <- if (k <= 4) ring_lab[k] else "park interior"
  r <- match_md(Tk, C7, COV_A)
  saveRDS(r$md, file.path(DIR_PROC, "matched", paste0("farring_ring", k, ".rds")))
  bt <- bal.tab(r$m, stats = c("m", "v"), un = TRUE, s.d.denom = "treated")$Balance; bt <- bt[rownames(bt) %in% COV_A, ]
  bal[[k]] <- data.table(run = "headline_vs_400_500m", ring = lab, covariate = rownames(bt), smd_before = bt$Diff.Un, smd_after = bt$Diff.Adj, vr_before = bt$V.Ratio.Un, vr_after = bt$V.Ratio.Adj)
  pr <- pairs_of(r$md); pr[, ring := lab]
  P[[k]] <- pr[, .(d = mean(lst_t - lst_c), n_scenes = .N), by = .(cell_t, park_id, block, ring)]
  logmsg(sprintf("Far-ring %-14s | %s | treated locations %d | trimmed %d | matched locations %d | mean diff %.2f C",
                 lab, r$method, nrow(Tk), r$trimmed, nrow(P[[k]]), mean(P[[k]]$d)))
}
fwrite(rbindlist(bal), file.path(DIR_OUT, "balance_table_headline.csv"))
L <- rbindlist(P)
## on-site (interior) headline estimate with block bootstrap
bootmean <- function(v, blk, B = NBOOT) { bl <- unique(blk); idx <- split(seq_along(v), blk)
  vapply(seq_len(B), function(b) { set.seed(SEED + b); mean(v[unlist(idx[sample(bl, length(bl), TRUE)], use.names = FALSE)]) }, 0) }
I <- L[ring == "park interior"]; bi <- bootmean(I$d, I$block)
logmsg(sprintf("Headline on-site (interior vs 400-500 m ring): %.2f C [%.2f, %.2f], %d locations", mean(I$d), quantile(bi, .025), quantile(bi, .975), nrow(I)))
fwrite(data.table(ring = "park interior", att = mean(I$d), se = sd(bi), ci_lo = quantile(bi, .025), ci_hi = quantile(bi, .975), locations = nrow(I)),
       file.path(DIR_OUT, "headline_onsite_vs_400_500m.csv"))

## 3. Per-park estimates by band
L <- L[ring != "park interior"]
L[, band := fifelse(ring == "0-50 m", "b0_50", fifelse(ring %in% c("50-100 m", "100-150 m"), "b50_150", "b150_200"))]
park_ci <- function(d, blk) {
  if (length(d) < 5) return(c(mean(d), NA, NA))
  bs <- if (uniqueN(blk) >= 3) bootmean(d, blk, 300) else vapply(seq_len(300), function(b) { set.seed(SEED + b); mean(sample(d, replace = TRUE)) }, 0)
  c(mean(d), quantile(bs, .025), quantile(bs, .975)) }
PP <- L[, { z <- park_ci(d, block); .(n_loc = .N, att = z[1], ci_lo = z[2], ci_hi = z[3]) }, by = .(park_id, band)]
W  <- dcast(PP, park_id ~ band, value.var = c("n_loc", "att", "ci_lo", "ci_hi"))
sig <- function(b) { x <- W[[paste0("ci_hi_", b)]]; !is.na(x) & x < 0 }
W[, cooling_distance_m := fifelse(sig("b0_50"), fifelse(sig("b50_150"), fifelse(sig("b150_200"), 200, 150), 50), 0)]
pa <- as.data.table(st_drop_geometry(parks))[, .(park_id = as.character(UNIT_ID), unit_name, area_ha, perim_m, pa_ratio, canopy_frac, water_ha, has_water, yr_est)]
W <- merge(pa, W, by = "park_id", all.x = FALSE)
fwrite(W, file.path(DIR_OUT, "per_park_att.csv"))
logmsg("Parks with any matched ring pixels: ", nrow(W), " | with >=5 locations at 50-150 m: ", sum(W$n_loc_b50_150 >= 5, na.rm = TRUE))
cat("\nCooling distance (contiguous significant bands), parks with >=5 locations in every band:\n")
print(W[n_loc_b0_50 >= 5 & n_loc_b50_150 >= 5 & n_loc_b150_200 >= 5, .N, by = cooling_distance_m][order(cooling_distance_m)])

## 4. Pooled 50-150 m cooling by park-area class
brks <- c(0.2, 0.5, 1, 2, 4, 8, 16, Inf)
L2 <- merge(L[band == "b50_150"], pa[, .(park_id, area_ha)], by = "park_id")
L2[, cls := cut(area_ha, brks, right = FALSE)]
AB <- L2[, { bs <- bootmean(d, block); .(parks = uniqueN(park_id), locations = .N, att = mean(d), ci_lo = quantile(bs, .025), ci_hi = quantile(bs, .975)) }, by = cls][order(cls)]
AB[, lower_ha := brks[as.integer(cls)]]; AB[, detect := ci_hi < 0]
print(AB)
thr <- NA; for (i in seq_len(nrow(AB))) if (all(AB$detect[i:nrow(AB)])) { thr <- AB$lower_ha[i]; break }
logmsg("Area threshold (cooling beyond 50 m detectable in this and every larger class): ", ifelse(is.na(thr), "none found", paste0(thr, " ha")))
fwrite(AB, file.path(DIR_OUT, "cooling_beyond_50m_by_park_area.csv"))

## 5. Meta-regression
M <- W[n_loc_b50_150 >= 5]
fit <- feols(att_b50_150 ~ log(area_ha) + pa_ratio + canopy_frac + has_water, data = M, weights = ~n_loc_b50_150, vcov = "hetero")
print(summary(fit)); cat("\nCorrelation log(area) vs P/A ratio:", round(cor(log(M$area_ha), M$pa_ratio), 2), "\n")
ct <- as.data.table(coeftable(fit), keep.rownames = "term"); fwrite(ct, file.path(DIR_OUT, "heterogeneity_meta_regression.csv"))
logmsg("Meta-regression parks: ", nobs(fit))
logmsg("STAGE_08_DONE")
