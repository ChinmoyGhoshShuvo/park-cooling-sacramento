# 04_park_screening.R
# Dissolves CPAD holdings into park units, screens them (area >= 0.2 ha,
# summer median NDVI > 0.3), records why each unit was kept or dropped, and
# cross-checks against OpenStreetMap park polygons. Runnable on its own.
source("code/_common.R")   # run from the repository folder
suppressPackageStartupMessages({ library(exactextractr); library(osmdata) })
set.seed(SEED)
ua   <- st_read(file.path(DIR_PROC, "study_area.gpkg"), quiet = TRUE)
nd   <- rast(file.path(DIR_PROC, "ndvi_median_summer.tif"))
cpad <- st_read(file.path(DIR_RAW, "cpad_holdings_sacramento.gpkg"), quiet = TRUE)
n_water <- sum(cpad$LAND_WATER == "Water", na.rm = TRUE)
cpad <- cpad[is.na(cpad$LAND_WATER) | cpad$LAND_WATER != "Water", ]
cpad$YR_EST[!is.na(cpad$YR_EST) & cpad$YR_EST <= 0] <- NA

## 1. dissolve holdings into park units
dt <- as.data.table(st_drop_geometry(cpad))
attrs <- dt[, .(unit_name = UNIT_NAME[1], agency_level = AGNCY_LEV[1], access = ACCESS_TYP[1],
                n_holdings = .N,
                yr_est = if (all(is.na(YR_EST))) NA_integer_ else as.integer(min(YR_EST, na.rm = TRUE))), by = UNIT_ID]
geo <- aggregate(st_geometry(cpad), by = list(UNIT_ID = cpad$UNIT_ID), FUN = st_union)
units <- merge(st_sf(geo), attrs, by = "UNIT_ID")
units <- st_make_valid(st_set_crs(units, st_crs(cpad)))
stopif(st_crs(units)$epsg != 32610, "park units lost their CRS")
units$area_ha   <- as.numeric(st_area(units)) / 1e4
units$ndvi_mean <- exact_extract(nd, units, "mean", progress = FALSE)
units$in_ua     <- lengths(st_intersects(units, ua)) > 0
units$status <- with(units, ifelse(!in_ua, "dropped_outside_study_area",
                     ifelse(area_ha < 0.2, "dropped_under_0.2ha",
                     ifelse(is.na(ndvi_mean), "dropped_no_ndvi",
                     ifelse(ndvi_mean <= 0.3, "dropped_ndvi_le_0.3", "kept")))))
st_write(units, file.path(DIR_PROC, "cpad_units_screened.gpkg"), delete_dsn = TRUE, quiet = TRUE)
parks <- units[units$status == "kept", ]
st_write(parks, file.path(DIR_PROC, "parks_kept.gpkg"), delete_dsn = TRUE, quiet = TRUE)

cat("\nCPAD water holdings removed before dissolving:", n_water, "\n")
cat("\nPark-unit screening:\n")
sm <- as.data.table(st_drop_geometry(units))[, .(units = .N, area_ha = round(sum(area_ha), 1)), by = status]
print(sm); fwrite(sm, file.path(DIR_OUT, "park_screening_summary.csv"))
cat("\nKept parks: area quantiles (ha):\n"); print(round(quantile(parks$area_ha, c(0, .25, .5, .75, .9, 1)), 2))
cat("\nKept parks by agency level:\n"); print(table(parks$agency_level))
logmsg("Kept parks with an establishment year: ", sum(!is.na(parks$yr_est)),
       " | established after 1975: ", sum(parks$yr_est > 1975, na.rm = TRUE))
stopif(nrow(parks) < 50, "fewer than 50 parks kept")

## 2. OpenStreetMap cross-check (report only; does not change the park set)
osm <- tryCatch({
  bb <- st_bbox(st_transform(ua, 4326))
  q  <- opq(bbox = as.numeric(bb), timeout = 300) |>
        add_osm_features(features = list("leisure" = c("park", "garden"), "landuse" = "recreation_ground"))
  o  <- osmdata_sf(q)
  g  <- rbind(o$osm_polygons[, "osm_id"], o$osm_multipolygons[, "osm_id"])
  st_make_valid(st_transform(g, CRS_WORK))
}, error = function(e) { logmsg("OSM cross-check FAILED: ", conditionMessage(e)); NULL })
if (!is.null(osm)) {
  osm <- osm[lengths(st_intersects(osm, ua)) > 0, ]
  osm$area_ha   <- as.numeric(st_area(osm)) / 1e4
  osm$ndvi_mean <- exact_extract(nd, osm, "mean", progress = FALSE)
  osm <- osm[osm$area_ha >= 0.2 & !is.na(osm$ndvi_mean) & osm$ndvi_mean > 0.3, ]
  cu  <- st_union(units)
  ov  <- as.numeric(st_area(st_intersection(osm, cu)))
  osm$overlap <- 0; idx <- unique(unlist(st_intersects(osm, cu)))
  osm$overlap <- vapply(seq_len(nrow(osm)), function(i) {
    x <- suppressWarnings(st_intersection(st_geometry(osm)[i], cu)); if (length(x) == 0) 0 else sum(as.numeric(st_area(x))) / (osm$area_ha[i] * 1e4) }, 0)
  only <- osm[osm$overlap < 0.10, ]
  st_write(only, file.path(DIR_PROC, "osm_parks_not_in_cpad.gpkg"), delete_dsn = TRUE, quiet = TRUE)
  logmsg("OSM green park polygons (>=0.2 ha, NDVI>0.3) in study area: ", nrow(osm), " (", round(sum(osm$area_ha)), " ha)")
  logmsg("  of which <10% covered by any CPAD land: ", nrow(only), " (", round(sum(only$area_ha)), " ha)")
  cat("\nLargest OSM-only parks (ha):\n"); print(head(sort(round(only$area_ha, 1), decreasing = TRUE), 15))
}
logmsg("STAGE_04_DONE")
