# 03_lst_ndvi.R
# Builds per-scene LST (deg C) and NDVI rasters on the native Landsat 30 m grid
# (EPSG:32610), cropped to the study area + 2 km. Runnable on its own.
source("code/_common.R")   # run from the repository folder
suppressPackageStartupMessages(library(rstac))
set.seed(SEED)
dir.create(file.path(DIR_PROC, "scenes"), showWarnings = FALSE)
ua  <- st_read(file.path(DIR_PROC, "study_area.gpkg"), quiet = TRUE)
sel <- fread(file.path(DIR_PROC, "landsat_scenes_selected.csv"))
aoi <- vect(st_buffer(st_union(ua), 2000))
uav <- vect(ua)

s  <- stac("https://planetarycomputer.microsoft.com/api/stac/v1")
it <- s |> stac_search(collections = "landsat-c2-l2", ids = sel$id, limit = 100) |> post_request() |> items_fetch()
stopif(length(it$features) != nrow(sel), paste("STAC returned", length(it$features), "of", nrow(sel), "scenes"))
ids <- vapply(it$features, `[[`, "", "id")
bitset <- function(x, b) bitwAnd(x, bitwShiftL(1L, b)) > 0
tmpl <- NULL; summ <- list()
for (k in seq_len(nrow(sel))) {
  if (k %% 8 == 1) signed <- items_sign(it, sign_fn = sign_planetary_computer())
  f  <- signed$features[[match(sel$id[k], ids)]]
  rd <- function(a) crop(rast(paste0("/vsicurl/", f$assets[[a]]$href)), ext(aoi), snap = "out")
  qa <- rd("qa_pixel"); st <- rd("lwir11"); b4 <- rd("red"); b5 <- rd("nir08")
  if (is.null(tmpl)) { tmpl <- rast(qa); values(tmpl) <- NA
                       writeRaster(tmpl, file.path(DIR_PROC, "grid_template.tif"), overwrite = TRUE) }
  for (r in list(qa, st, b4, b5)) stopif(!compareGeom(r, tmpl, stopOnError = FALSE), paste(sel$id[k], "grid differs from template"))
  q   <- values(qa, mat = FALSE)
  bad <- bitset(q, 0) | bitset(q, 1) | bitset(q, 2) | bitset(q, 3) | bitset(q, 4) | bitset(q, 7)
  stv <- values(st, mat = FALSE); stv[stv == 0] <- NA
  lst <- stv * 0.00341802 + 149.0 - 273.15; lst[bad] <- NA
  red <- values(b4, mat = FALSE) * 0.0000275 - 0.2; nir <- values(b5, mat = FALSE) * 0.0000275 - 0.2
  ndvi <- (nir - red) / (nir + red); ndvi[bad | !is.finite(ndvi) | abs(ndvi) > 1] <- NA   # reflectance can go slightly negative over dark pixels
  L <- setValues(rast(tmpl), lst); N <- setValues(rast(tmpl), ndvi)
  d <- format(as.Date(sel$date[k]), "%Y%m%d")
  writeRaster(L, file.path(DIR_PROC, "scenes", paste0("LST_", d, ".tif")), overwrite = TRUE, datatype = "FLT4S")
  writeRaster(N, file.path(DIR_PROC, "scenes", paste0("NDVI_", d, ".tif")), overwrite = TRUE, datatype = "FLT4S")
  Lu <- values(mask(L, uav), mat = FALSE); Nu <- values(mask(N, uav), mat = FALSE)
  inua <- !is.na(values(mask(setValues(rast(tmpl), 1), uav), mat = FALSE))
  summ[[k]] <- data.table(date = d, platform = sel$platform[k],
    lst_p01 = quantile(Lu, .01, na.rm = TRUE), lst_med = median(Lu, na.rm = TRUE), lst_p99 = quantile(Lu, .99, na.rm = TRUE),
    lst_min = min(Lu, na.rm = TRUE), lst_max = max(Lu, na.rm = TRUE),
    ndvi_med = median(Nu, na.rm = TRUE), ndvi_min = min(Nu, na.rm = TRUE), ndvi_max = max(Nu, na.rm = TRUE),
    na_frac_ua = mean(is.na(Lu[inua])))
  logmsg("scene ", k, "/", nrow(sel), " ", d, " done")
}
S <- rbindlist(summ); fwrite(S, file.path(DIR_PROC, "scene_lst_ndvi_summary.csv"))
cat("\nTemplate grid:\n"); print(tmpl)
cat("\nPer-scene summary (study area only):\n"); print(S[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
stopif(any(S$lst_med < 25 | S$lst_med > 55), "median summer LST outside 25-55 C")
stopif(any(S$lst_min < 5 | S$lst_max > 80), "LST extremes implausible")
stopif(any(S$ndvi_min < -1.01 | S$ndvi_max > 1.01), "NDVI outside [-1, 1]")
stopif(any(S$na_frac_ua > 0.15), "more than 15% of study area masked in a scene")

## Summer-mean NDVI composite (used for park screening)
nd <- rast(list.files(file.path(DIR_PROC, "scenes"), "^NDVI_.*tif$", full.names = TRUE))
mn <- app(nd, median, na.rm = TRUE); writeRaster(mn, file.path(DIR_PROC, "ndvi_median_summer.tif"), overwrite = TRUE)
logmsg("NDVI summer-median composite written: ", nlyr(nd), " scenes")
logmsg("STAGE_03_DONE")
