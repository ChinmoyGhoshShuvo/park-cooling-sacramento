# 03b_ndvi_clean_composite.R
# Post-processes the stage 03 NDVI rasters: blanks pixels with |NDVI| > 1
# (dark pixels where scaled reflectance goes negative), validates all 20 scenes,
# and builds the summer median NDVI composite. Runnable on its own.
source("code/_common.R")   # run from the repository folder
set.seed(SEED)
sel  <- fread(file.path(DIR_PROC, "landsat_scenes_selected.csv"))
d    <- format(as.Date(sel$date), "%Y%m%d")
tmpl <- rast(file.path(DIR_PROC, "grid_template.tif"))
ua   <- vect(st_read(file.path(DIR_PROC, "study_area.gpkg"), quiet = TRUE))
stopif(crs(tmpl, describe = TRUE)$code != "32610" || any(res(tmpl) != 30), "template grid not 30 m EPSG:32610")
out <- list()
for (k in seq_along(d)) {
  fl <- file.path(DIR_PROC, "scenes", paste0("LST_", d[k], ".tif")); fn <- file.path(DIR_PROC, "scenes", paste0("NDVI_", d[k], ".tif"))
  stopif(!file.exists(fl) || !file.exists(fn), paste("missing scene files for", d[k]))
  L <- rast(fl); N <- rast(fn)
  stopif(!compareGeom(L, tmpl, stopOnError = FALSE) || !compareGeom(N, tmpl, stopOnError = FALSE), paste(d[k], "grid mismatch"))
  v <- values(N, mat = FALSE); badn <- !is.na(v) & abs(v) > 1
  if (any(badn)) { v[badn] <- NA; N2 <- setValues(rast(tmpl), v)
    writeRaster(N2, paste0(fn, ".tmp.tif"), overwrite = TRUE, datatype = "FLT4S"); rm(N); gc()
    file.remove(fn); file.rename(paste0(fn, ".tmp.tif"), fn) }
  out[[k]] <- data.table(date = d[k], ndvi_px_blanked = sum(badn), share_of_valid = signif(sum(badn) / sum(!is.na(v) | badn), 3))
}
O <- rbindlist(out); fwrite(O, file.path(DIR_PROC, "ndvi_outliers_blanked.csv")); print(O)
nd <- rast(file.path(DIR_PROC, "scenes", paste0("NDVI_", d, ".tif")))
rg <- minmax(nd); stopif(any(rg < -1 | rg > 1, na.rm = TRUE), "NDVI still outside [-1, 1]")
med <- app(nd, median, na.rm = TRUE)
writeRaster(med, file.path(DIR_PROC, "ndvi_median_summer.tif"), overwrite = TRUE, datatype = "FLT4S")
mv <- values(mask(med, ua), mat = FALSE)
logmsg("Median NDVI composite in study area: p05/p50/p95 = ", paste(round(quantile(mv, c(.05, .5, .95), na.rm = TRUE), 3), collapse = "/"),
       " | NA ", round(100 * mean(is.na(mv[!is.na(values(mask(setValues(rast(tmpl), 1), ua), mat = FALSE))])), 2), "%")
logmsg("Total NDVI pixels blanked across 20 scenes: ", sum(O$ndvi_px_blanked))
logmsg("STAGE_03b_DONE")
