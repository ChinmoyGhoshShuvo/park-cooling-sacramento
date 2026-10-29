# 05c_covariate_ghsl1975.R
# GHSL GHS-BUILT-S 1975 (100 m, Mollweide) -> built fraction on the Landsat grid.
# Used ONLY as the negative-control outcome. Runnable on its own.
source("code/_common.R")   # run from the repository folder
set.seed(SEED)
tmpl <- rast(file.path(DIR_PROC, "grid_template.tif"))
ua   <- st_read(file.path(DIR_PROC, "study_area.gpkg"), quiet = TRUE); uav <- vect(ua)
inua <- !is.na(values(mask(setValues(rast(tmpl), 1), uav), mat = FALSE))
dir.create(file.path(DIR_RAW, "ghsl"), showWarnings = FALSE)
bbm <- st_bbox(st_transform(st_as_sfc(st_bbox(tmpl)), "ESRI:54009"))
cols <- unique(floor((c(bbm["xmin"], bbm["xmax"]) + 18041000) / 1e6) + 1)
rows <- unique(floor((9e6 - c(bbm["ymax"], bbm["ymin"])) / 1e6) + 1)
base <- "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A/GHS_BUILT_S_E1975_GLOBE_R2023A_54009_100/V1-0/tiles/"
tifs <- c()
for (R in rows) for (C in cols) {
  nm <- sprintf("GHS_BUILT_S_E1975_GLOBE_R2023A_54009_100_V1_0_R%d_C%d", R, C)
  z  <- file.path(DIR_RAW, "ghsl", paste0(nm, ".zip"))
  if (!file.exists(z)) download.file(paste0(base, nm, ".zip"), z, mode = "wb")
  unzip(z, exdir = file.path(DIR_RAW, "ghsl"), overwrite = TRUE)
  tifs <- c(tifs, file.path(DIR_RAW, "ghsl", paste0(nm, ".tif")))
  logmsg("GHSL tile ", nm)
}
g <- if (length(tifs) > 1) vrt(tifs) else rast(tifs)
stopif(!relate(ext(g), ext(as.numeric(bbm[c("xmin","xmax","ymin","ymax")])), "contains"), "GHSL tiles do not cover the study window")
g <- crop(g, ext(as.numeric(bbm[c("xmin","xmax","ymin","ymax")])) + 500)
g[g > 10000] <- NA
b75 <- project(g, tmpl, method = "bilinear") / 10000; names(b75) <- "built1975"
v <- values(b75, mat = FALSE)[inua]; q <- quantile(v, c(0, .5, .9, .99, 1), na.rm = TRUE)
logmsg(sprintf("built1975  NA %.3f%% | min %.3f med %.3f p90 %.3f p99 %.3f max %.3f | share of pixels with any 1975 built-up %.1f%%",
               100 * mean(is.na(v)), q[1], q[2], q[3], q[4], q[5], 100 * mean(v > 0.01, na.rm = TRUE)))
stopif(!compareGeom(b75, tmpl, stopOnError = FALSE), "grid mismatch")
stopif(mean(is.na(v)) > 0.001, "GHSL has gaps in the study area")
stopif(q[1] < 0 | q[5] > 1.0001, "built fraction outside [0, 1]")
writeRaster(b75, file.path(DIR_PROC, "cov", "built1975.tif"), overwrite = TRUE, datatype = "FLT4S")
logmsg("STAGE_05c_DONE")
