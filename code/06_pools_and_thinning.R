# 06_pools_and_thinning.R
# Defines treated (distance rings from park edge), control (>800 m from any
# CPAD/OSM green space and from permanent water) and park-interior pixels for
# the MAIN (urban parks) and SECONDARY (all green CPAD) definitions; thins all
# pools on one 120 m lattice; attaches covariates and per-scene LST. Runnable on its own.
source("code/_common.R")   # run from the repository folder
set.seed(SEED)
tmpl <- rast(file.path(DIR_PROC, "grid_template.tif"))
ua   <- vect(st_read(file.path(DIR_PROC, "study_area.gpkg"), quiet = TRUE))
inua <- !is.na(values(mask(setValues(rast(tmpl), 1), ua), mat = FALSE))
pm   <- vect(file.path(DIR_PROC, "parks_main_urban.gpkg"))
pa   <- vect(file.path(DIR_PROC, "parks_secondary_all_green.gpkg"))
ex   <- vect(file.path(DIR_PROC, "park_exclusion_all.gpkg"))
perm <- values(rast(file.path(DIR_PROC, "cov", "perm_water.tif")), mat = FALSE) == 1
dwat <- values(rast(file.path(DIR_PROC, "cov", "dist_water_m.tif")), mat = FALSE)
RINGS <- c(0, 50, 100, 150, 200, 300, 400, 500)

## distance from each cell centre to the nearest park polygon (0 inside)
t0 <- Sys.time()
dist_to <- function(v, nm) {
  inside <- rasterize(v, rast(tmpl), field = 1, background = 0, touches = FALSE)
  bnd <- as.lines(v)
  d <- distance(rast(tmpl), bnd)                       # distance to nearest boundary line
  d <- ifel(inside == 1, 0, d); names(d) <- nm
  writeRaster(d, file.path(DIR_PROC, "cov", paste0(nm, ".tif")), overwrite = TRUE, datatype = "FLT4S")
  list(d = values(d, mat = FALSE), inside = values(inside, mat = FALSE) == 1, insr = inside)
}
DM <- dist_to(pm, "dist_park_main"); DA <- dist_to(pa, "dist_park_all"); DX <- dist_to(ex, "dist_any_green")
logmsg("Distance rasters built in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
for (z in list(DM$d, DA$d, DX$d)) stopif(any(is.na(z[inua])) || min(z[inua]) < 0, "bad park-distance values")

## pure park interior: cell and its 8 neighbours all inside a park
interior <- function(r) values(focal(r, w = 3, fun = "min", na.policy = "omit", fillvalue = 0), mat = FALSE) == 1
IM <- interior(DM$insr); IA <- interior(DA$insr)

land <- inua & !perm
ring <- function(d) as.integer(cut(d, RINGS, right = TRUE, include.lowest = FALSE))
pool <- data.table(cell = seq_along(inua), land, d_main = DM$d, d_all = DA$d, d_any = DX$d, dwat,
                   in_any = DX$inside, int_main = IM, int_all = IA)
pool[, ring_main := fifelse(land & !in_any & d_main > 0 & d_main <= 500, ring(d_main), NA_integer_)]
pool[, ring_all  := fifelse(land & !in_any & d_all  > 0 & d_all  <= 500, ring(d_all),  NA_integer_)]
pool[, control   := land & !in_any & d_any > 800 & dwat > 800]
pool[, interior_main := inua & !perm & int_main]
pool[, interior_all  := inua & !perm & int_all]
pool <- pool[!is.na(ring_main) | !is.na(ring_all) | control | interior_main | interior_all]
logmsg("Unthinned pools | treated main ", sum(!is.na(pool$ring_main)), " | treated all ", sum(!is.na(pool$ring_all)),
       " | control ", sum(pool$control), " | interior main ", sum(pool$interior_main))

## 120 m thinning: one lattice (every 4th row and column) with a seeded random origin
xy <- xyFromCell(tmpl, pool$cell); rc <- rowColFromCell(tmpl, pool$cell)
off <- sample(0:3, 2); logmsg("Thinning lattice origin offset (row, col): ", paste(off, collapse = ", "))
pool[, `:=`(x = xy[, 1], y = xy[, 2], keep = ((rc[, 1] - off[1]) %% 4 == 0) & ((rc[, 2] - off[2]) %% 4 == 0))]
pool <- pool[keep == TRUE][, keep := NULL]
pool[, block := paste0(floor(x / 1000), "_", floor(y / 1000))]
md <- min(dist(as.matrix(pool[sample(.N, min(.N, 3000)), .(x, y)])))
logmsg("Min spacing in a 3,000-pixel sample after thinning: ", md, " m"); stopif(md < 120, "thinning failed")

## nearest park id for treated pixels
pts <- st_as_sf(pool[, .(x, y)], coords = c("x", "y"), crs = CRS_WORK)
pms <- st_as_sf(pm); pas <- st_as_sf(pa)
pool[, park_main := NA_character_]; pool[, park_all := NA_character_]
i <- which(!is.na(pool$ring_main) | pool$interior_main); pool$park_main[i] <- as.character(pms$UNIT_ID[st_nearest_feature(pts[i, ], pms)])
i <- which(!is.na(pool$ring_all)  | pool$interior_all);  pool$park_all[i]  <- as.character(pas$UNIT_ID[st_nearest_feature(pts[i, ], pas)])

## covariates
cv <- c("elev_m", "slope_deg", "northness", "eastness", "dist_water_m", "dist_cbd_m", "dist_coast_m", "bldg_dens150", "built1975")
for (v in cv) {
  f <- file.path(DIR_PROC, "cov", paste0(v, ".tif")); stopif(!file.exists(f), paste("missing covariate", v))
  pool[[v]] <- values(rast(f), mat = FALSE)[pool$cell]
}
stopif(any(is.na(pool[, ..cv])), "NA in covariates for pooled pixels")

## per-scene LST (long format)
sel <- fread(file.path(DIR_PROC, "landsat_scenes_selected.csv")); d <- format(as.Date(sel$date), "%Y%m%d")
long <- rbindlist(lapply(d, function(s) data.table(cell = pool$cell, scene = s,
          lst = values(rast(file.path(DIR_PROC, "scenes", paste0("LST_", s, ".tif"))), mat = FALSE)[pool$cell])))
long <- long[!is.na(lst)]
saveRDS(pool, file.path(DIR_PROC, "pool_locations.rds")); saveRDS(long, file.path(DIR_PROC, "pool_lst_long.rds"))

cat("\nThinned locations per ring (main | all green):\n")
print(merge(pool[!is.na(ring_main), .(main = .N), by = .(ring = ring_main)], pool[!is.na(ring_all), .(all_green = .N), by = .(ring = ring_all)], by = "ring")[
  , ring := paste0(RINGS[ring], "-", RINGS[ring + 1], " m")][])
logmsg("Thinned | treated main ", sum(!is.na(pool$ring_main)), " | treated all ", sum(!is.na(pool$ring_all)),
       " | control ", sum(pool$control), " | interior main ", sum(pool$interior_main), " | interior all ", sum(pool$interior_all),
       " | blocks ", uniqueN(pool$block))
logmsg("Pixel-scene LST observations: ", nrow(long), " (", round(100 * nrow(long) / (nrow(pool) * length(d)), 2), "% of possible)")
stopif(sum(pool$control) < 1000, "fewer than 1,000 control locations")
logmsg("STAGE_06_DONE")
