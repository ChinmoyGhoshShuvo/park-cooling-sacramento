# 05b_covariate_buildings.R
# Present-day built density (SECONDARY specification only): share of land under
# building footprints within a 150 m radius, from Microsoft US Building
# Footprints (California). Runnable on its own.
source("code/_common.R")   # run from the repository folder
set.seed(SEED)
tmpl <- rast(file.path(DIR_PROC, "grid_template.tif"))
ua   <- st_read(file.path(DIR_PROC, "study_area.gpkg"), quiet = TRUE); uav <- vect(ua)
inua <- !is.na(values(mask(setValues(rast(tmpl), 1), uav), mat = FALSE))
dir.create(file.path(DIR_PROC, "cov"), showWarnings = FALSE)
zipf <- file.path(DIR_RAW, "California.geojson.zip")
if (!file.exists(zipf)) download.file("https://minedbuildings.z5.web.core.windows.net/legacy/usbuildings-v2/California.geojson.zip", zipf, mode = "wb")
inner <- unzip(zipf, list = TRUE)$Name; inner <- inner[grepl("geojson$", inner)][1]
out <- file.path(DIR_RAW, "ms_buildings_sacramento.gpkg")
bb <- as.numeric(st_bbox(st_transform(st_as_sfc(st_bbox(tmpl)), 4326)))
if (!file.exists(out)) {
  logmsg("Extracting buildings for the study window from ", inner, " (streams the whole state file)")
  sf::gdal_utils("vectortranslate", paste0("/vsizip/", zipf, "/", inner), out,
                 options = c("-spat", bb, "-t_srs", "EPSG:32610", "-nln", "buildings", "-f", "GPKG", "-overwrite"))
}
b <- vect(out); logmsg("Buildings in window: ", nrow(b))
stopif(nrow(b) < 100000, "suspiciously few buildings extracted")
cov <- rasterize(b, rast(tmpl), cover = TRUE, background = 0)
cov <- clamp(cov, 0, 1)
w <- focalMat(cov, 150, "circle")
dens <- focal(cov, w, fun = "sum", na.rm = TRUE); names(dens) <- "bldg_dens150"
v <- values(dens, mat = FALSE)[inua]; q <- quantile(v, c(0, .5, .9, .99, 1), na.rm = TRUE)
logmsg(sprintf("bldg_dens150  NA %.3f%% | min %.3f med %.3f p90 %.3f p99 %.3f max %.3f", 100 * mean(is.na(v)), q[1], q[2], q[3], q[4], q[5]))
stopif(!compareGeom(dens, tmpl, stopOnError = FALSE), "grid mismatch")
stopif(q[1] < 0 | q[5] > 1.0001, "built fraction outside [0, 1]")
stopif(q[2] < 0.03 | q[2] > 0.5, "median built fraction implausible for a US urban area")
writeRaster(dens, file.path(DIR_PROC, "cov", "bldg_dens150.tif"), overwrite = TRUE, datatype = "FLT4S")
logmsg("STAGE_05b_DONE")
