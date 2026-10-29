# 05a_covariates_terrain_water.R
# Exogenous covariates on the Landsat grid: elevation, slope, aspect (as
# northness/eastness), distance to permanent water (JRC GSW occurrence > 50%),
# distance to CBD, distance to coastline (Natural Earth 1:10m). Runnable on its own.
source("code/_common.R")   # run from the repository folder
suppressPackageStartupMessages(library(rstac))
set.seed(SEED)
tmpl <- rast(file.path(DIR_PROC, "grid_template.tif"))
ua   <- st_read(file.path(DIR_PROC, "study_area.gpkg"), quiet = TRUE); uav <- vect(ua)
inua <- !is.na(values(mask(setValues(rast(tmpl), 1), uav), mat = FALSE))
dir.create(file.path(DIR_PROC, "cov"), showWarnings = FALSE)
bb   <- as.numeric(st_bbox(st_transform(st_as_sfc(st_bbox(tmpl)), 4326)))
s    <- stac("https://planetarycomputer.microsoft.com/api/stac/v1")
wopt <- list(datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
check <- function(r, name, lo, hi, max_na = 0.001) {
  stopif(!compareGeom(r, tmpl, stopOnError = FALSE), paste(name, "grid mismatch"))
  v <- values(r, mat = FALSE)[inua]; na <- mean(is.na(v)); q <- quantile(v, c(0, .01, .5, .99, 1), na.rm = TRUE)
  logmsg(sprintf("%-14s NA %.3f%% | min %.2f p01 %.2f med %.2f p99 %.2f max %.2f", name, 100 * na, q[1], q[2], q[3], q[4], q[5]))
  stopif(na > max_na, paste(name, "NA fraction too high"))
  stopif(q[1] < lo | q[5] > hi, paste(name, "values outside plausible range", lo, "to", hi))
  writeRaster(r, file.path(DIR_PROC, "cov", paste0(name, ".tif")), overwrite = TRUE, wopt = wopt)
}

## Elevation (Copernicus GLO-30), bilinear onto the Landsat grid
dem_it <- s |> stac_search(collections = "cop-dem-glo-30", bbox = bb, limit = 20) |> post_request() |> items_sign(sign_fn = sign_planetary_computer())
logmsg("DEM tiles: ", length(dem_it$features))
dem_src <- vrt(paste0("/vsicurl/", vapply(dem_it$features, function(f) f$assets$data$href, "")))
dem <- project(dem_src, tmpl, method = "bilinear"); names(dem) <- "elev"
check(dem, "elev_m", -10, 700)
tr <- terrain(dem, v = c("slope", "aspect"), unit = "degrees")
check(tr$slope, "slope_deg", 0, 75)
asp <- tr$aspect
north <- cos(asp * pi / 180); east <- sin(asp * pi / 180)
north[tr$slope < 0.5] <- 0; east[tr$slope < 0.5] <- 0      # flat ground has no meaningful aspect
north[is.na(north) & !is.na(dem)] <- 0; east[is.na(east) & !is.na(dem)] <- 0
check(north, "northness", -1, 1); check(east, "eastness", -1, 1)

## Permanent water (JRC GSW occurrence > 50%)
gsw_it <- s |> stac_search(collections = "jrc-gsw", bbox = bb, limit = 20) |> post_request() |> items_sign(sign_fn = sign_planetary_computer())
logmsg("JRC GSW tiles: ", length(gsw_it$features))
occ <- vrt(paste0("/vsicurl/", vapply(gsw_it$features, function(f) f$assets$occurrence$href, "")))
occ <- crop(occ, ext(bb[c(1, 3, 2, 4)]) + 0.02)
wat <- classify(occ, rbind(c(-Inf, 50, 0), c(50, 100, 1), c(100, Inf, NA)), include.lowest = TRUE, right = TRUE)
wfrac <- project(wat, tmpl, method = "average")
perm <- wfrac >= 0.5; names(perm) <- "perm_water"
writeRaster(perm, file.path(DIR_PROC, "cov", "perm_water.tif"), overwrite = TRUE, datatype = "INT1U")
pw <- values(perm, mat = FALSE)[inua]
logmsg("Permanent water: ", round(100 * mean(pw == 1, na.rm = TRUE), 2), "% of study area pixels")
stopif(mean(pw == 1, na.rm = TRUE) < 0.002 | mean(pw == 1, na.rm = TRUE) > 0.15, "permanent-water share implausible")
dw <- distance(classify(perm, cbind(0, NA))); names(dw) <- "dist_water"
check(dw, "dist_water_m", 0, 20000)

## Distance to CBD
cbd <- vect(st_read(file.path(DIR_PROC, "cbd.gpkg"), quiet = TRUE))
dc <- distance(rast(tmpl), cbd); names(dc) <- "dist_cbd"
check(dc, "dist_cbd_m", 0, 60000)

## Distance to coastline (Natural Earth 1:10m), computed at 300 m then interpolated
cz <- file.path(DIR_RAW, "ne_10m_coastline.geojson")
if (!file.exists(cz)) download.file("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_coastline.geojson", cz, mode = "wb")
coast <- st_read(cz, quiet = TRUE)
coast <- suppressWarnings(st_crop(coast, xmin = -126, xmax = -118, ymin = 34, ymax = 42))
coast <- st_transform(coast, CRS_WORK)
coarse <- aggregate(rast(tmpl), 10)
dco <- resample(distance(coarse, vect(coast)), tmpl, method = "bilinear"); names(dco) <- "dist_coast"
check(dco, "dist_coast_m", 30000, 250000)
cbd_sf <- st_read(file.path(DIR_PROC, "cbd.gpkg"), quiet = TRUE)
np <- st_cast(st_nearest_points(cbd_sf, st_union(coast)), "POINT")[2]
logmsg("Nearest Natural Earth coastline to the CBD: ", paste(round(st_coordinates(st_transform(np, 4326)), 3), collapse = ", "),
       " (lon, lat) at ", round(as.numeric(st_distance(cbd_sf, np)) / 1000, 1), " km")
logmsg("STAGE_05a_DONE")
