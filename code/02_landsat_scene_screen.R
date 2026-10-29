# 02_landsat_scene_screen.R
# Searches Landsat 8/9 C2 L2 summer scenes over the study area, screens each one
# for clear fraction, coverage and smoke (high aerosol), and applies the
# pre-registered selection rule. Runnable on its own.
source("code/_common.R")   # run from the repository folder
suppressPackageStartupMessages(library(rstac))
set.seed(SEED)
ua  <- st_read(file.path(DIR_PROC, "study_area.gpkg"), quiet = TRUE)
bb  <- as.numeric(st_bbox(st_transform(ua, 4326)))
s   <- stac("https://planetarycomputer.microsoft.com/api/stac/v1")
it  <- s |> stac_search(collections = "landsat-c2-l2", bbox = bb,
                        datetime = "2016-06-01T00:00:00Z/2025-09-30T23:59:59Z", limit = 250) |>
       ext_query(`eo:cloud_cover` < 30) |> post_request() |> items_fetch()
p   <- lapply(it$features, `[[`, "properties")
cand <- data.table(id = vapply(it$features, `[[`, "", "id"),
                   date = as.Date(substr(vapply(p, `[[`, "", "datetime"), 1, 10)),
                   platform = vapply(p, `[[`, "", "platform"),
                   path = vapply(p, function(x) as.integer(x[["landsat:wrs_path"]]), 1L),
                   row  = vapply(p, function(x) as.integer(x[["landsat:wrs_row"]]), 1L),
                   tier = vapply(p, function(x) x[["landsat:collection_category"]], ""),
                   epsg = vapply(p, function(x) as.integer(x[["proj:epsg"]] %||% NA), 1L),
                   cc   = vapply(p, `[[`, 1, "eo:cloud_cover"))
cand$idx <- seq_len(nrow(cand))
cand <- cand[platform %in% c("landsat-8", "landsat-9") & tier == "T1" &
             as.integer(format(date, "%m")) %in% 6:9]
logmsg("Candidates after platform/tier/month filter: ", nrow(cand))
print(cand[, .N, by = .(path, row, epsg)])

uav <- vect(st_transform(ua, CRS_WORK))
bitset <- function(x, b) bitwAnd(x, bitwShiftL(1L, b)) > 0
## Footprint pre-check: the tile outline must contain the whole study area
ua_ll <- st_union(st_transform(ua, 4326))
fp_ok <- vapply(cand$idx, function(i) {
  g <- st_as_sfc(list(it$features[[i]]$geometry) |> lapply(function(x) {
         x$coordinates <- lapply(x$coordinates, function(r) do.call(rbind, lapply(r, unlist)))
         if (x$type == "Polygon") st_polygon(x$coordinates) else st_multipolygon(list(x$coordinates)) }), crs = 4326)
  as.logical(st_contains(g, ua_ll, sparse = FALSE)[1, 1])
}, TRUE)
logmsg("Tile footprint contains the whole study area: ", sum(fp_ok), " of ", nrow(cand))
print(cand[, .(contain = sum(fp_ok[.I])), by = .(path, row)])
cand <- cand[fp_ok]

signed <- items_sign(it, sign_fn = sign_planetary_computer())
res <- vector("list", nrow(cand))
for (k in seq_len(nrow(cand))) {
  i  <- cand$idx[k]
  if (k %% 40 == 0) signed <- items_sign(it, sign_fn = sign_planetary_computer())   # refresh access tokens
  f  <- signed$features[[i]]
  out <- tryCatch({
    qa <- rast(paste0("/vsicurl/", f$assets$qa_pixel$href))
    if (!same.crs(qa, CRS_WORK)) stop("scene CRS ", crs(qa, describe = TRUE)$code)
    qa <- mask(crop(qa, uav), uav); v <- values(qa, mat = FALSE); v <- v[!is.na(v)]
    n_all <- length(v); fill <- bitset(v, 0)
    bad   <- fill | bitset(v, 1) | bitset(v, 2) | bitset(v, 3) | bitset(v, 4)
    ae <- rast(paste0("/vsicurl/", f$assets$qa_aerosol$href))
    ae <- mask(crop(ae, uav), uav); a <- values(ae, mat = FALSE); a <- a[!is.na(a) & !bitset(a, 0)]
    hi <- bitwAnd(bitwShiftR(a, 6L), 3L) == 3L
    data.table(coverage = 1 - mean(fill), clear = 1 - mean(bad), high_aerosol = mean(hi), n_px = n_all)
  }, error = function(e) data.table(coverage = NA_real_, clear = NA_real_, high_aerosol = NA_real_, n_px = NA_integer_, err = conditionMessage(e)))
  res[[k]] <- cbind(cand[k], out)
  if (k %% 10 == 0) logmsg("screened ", k, " / ", nrow(cand))
}
scr <- rbindlist(res, fill = TRUE)
fwrite(scr, file.path(DIR_PROC, "landsat_scene_screening.csv"))
if ("err" %in% names(scr)) { cat("\nScenes that failed to read:\n"); print(scr[!is.na(err), .(id, err)]) }

## Pre-registered selection rule
ok <- scr[!is.na(clear) & coverage >= 0.99 & clear >= 0.85 & high_aerosol <= 0.10]
ok[, yr := as.integer(format(date, "%Y"))]
ok[, dist_jul31 := abs(as.numeric(date - as.Date(paste0(yr, "-07-31"))))]
setorder(ok, yr, dist_jul31, -clear)
sel <- ok[, head(.SD, 2), by = yr]
fwrite(sel, file.path(DIR_PROC, "landsat_scenes_selected.csv"))

logmsg("Screened: ", nrow(scr), " | pass coverage>=99%: ", sum(scr$coverage >= .99, na.rm = TRUE),
       " | + clear>=85%: ", sum(scr$coverage >= .99 & scr$clear >= .85, na.rm = TRUE),
       " | + aerosol<=10%: ", nrow(ok))
cat("\nRejected for smoke (pass cloud + coverage, fail aerosol):\n")
print(scr[coverage >= .99 & clear >= .85 & high_aerosol > .10, .(id, date, clear = round(clear, 3), high_aerosol = round(high_aerosol, 3))])
cat("\nSELECTED SCENES:\n"); print(sel[, .(id, date, platform, path, row, clear = round(clear, 3), high_aerosol = round(high_aerosol, 3))])
stopif(nrow(sel) < 6, "fewer than 6 qualifying scenes")
stopif(length(unique(sel$yr)) < 4, "selected scenes span fewer than 4 years")
stopif(length(unique(paste(sel$path, sel$row))) > 1, "selected scenes come from more than one WRS path/row")
logmsg("STAGE_02_DONE")
