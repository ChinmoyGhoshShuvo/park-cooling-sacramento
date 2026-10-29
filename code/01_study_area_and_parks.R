# 01_study_area_and_parks.R
# Builds the Sacramento study area (2020 Census urban area), the CBD point,
# and downloads CPAD park/open-space holdings for the area. Runnable on its own.
source("code/_common.R")   # run from the repository folder
set.seed(SEED)

## 1. Study area: 2020 Census urban area "Sacramento, CA"
zip <- file.path(DIR_RAW, "tl_2023_us_uac20.zip")
if (!file.exists(zip)) download.file("https://www2.census.gov/geo/tiger/TIGER2023/UAC/tl_2023_us_uac20.zip", zip, mode = "wb")
ua_all <- st_read(paste0("/vsizip/", zip), quiet = TRUE)
logmsg("Urban-area fields: ", paste(names(ua_all), collapse = ", "))
namecol <- intersect(c("NAME20", "NAMELSAD20", "NAME"), names(ua_all))[1]
ua <- ua_all[grepl("^Sacramento, CA", ua_all[[namecol]]), ]
stopif(nrow(ua) != 1, paste("expected 1 Sacramento UA, found", nrow(ua)))
ua <- st_make_valid(st_transform(ua, CRS_WORK))
ua_km2 <- as.numeric(st_area(ua)) / 1e6
logmsg("Study area: ", ua[[namecol]], " | CRS ", st_crs(ua)$epsg, " | area ", round(ua_km2, 1), " km2")
logmsg("Study area bbox (m): ", paste(round(st_bbox(ua)), collapse = ", "))
stopif(st_crs(ua)$epsg != 32610, "study area not in EPSG:32610")
stopif(ua_km2 < 500 | ua_km2 > 2500, "study area size implausible")
st_write(ua, file.path(DIR_PROC, "study_area.gpkg"), delete_dsn = TRUE, quiet = TRUE)

## 2. CBD point
cbd <- st_transform(st_sfc(st_point(CBD_LONLAT), crs = 4326), CRS_WORK)
stopif(!as.logical(st_intersects(cbd, ua, sparse = FALSE)), "CBD point falls outside study area")
st_write(st_sf(name = "California State Capitol", geometry = cbd), file.path(DIR_PROC, "cbd.gpkg"), delete_dsn = TRUE, quiet = TRUE)

## 3. CPAD holdings (official CNRA map service), study area bbox + 2 km
bb <- st_bbox(st_transform(st_buffer(ua, 2000), 4326))
base <- "https://gis.cnra.ca.gov/arcgis/rest/services/Boundaries/CPAD_AgencyLevel/MapServer/0/query"
off <- 0; parts <- list(); page <- 500
repeat {
  q <- c(where = "1=1", geometry = paste(round(bb[c("xmin","ymin","xmax","ymax")], 6), collapse = ","),
         geometryType = "esriGeometryEnvelope", inSR = "4326", spatialRel = "esriSpatialRelIntersects",
         outFields = "*", outSR = "4326", f = "geojson", orderByFields = "OBJECTID",
         resultOffset = off, resultRecordCount = page)
  url <- paste0(base, "?", paste(names(q), vapply(q, URLencode, "", reserved = TRUE), sep = "=", collapse = "&"))
  tf <- tempfile(fileext = ".geojson"); download.file(url, tf, quiet = TRUE, mode = "wb")
  x <- st_read(tf, quiet = TRUE)
  if (nrow(x) == 0) break
  parts[[length(parts) + 1]] <- x; off <- off + nrow(x)
  logmsg("CPAD page fetched: ", nrow(x), " (total ", off, ")")
  if (nrow(x) < page) break
}
cpad <- do.call(rbind, parts)
stopif(anyDuplicated(cpad$OBJECTID) > 0, "duplicate CPAD records from paging")
cpad <- st_make_valid(st_transform(cpad, CRS_WORK))
cpad <- cpad[lengths(st_intersects(cpad, st_buffer(ua, 2000))) > 0, ]
cpad$area_ha <- as.numeric(st_area(cpad)) / 1e4
st_write(cpad, file.path(DIR_RAW, "cpad_holdings_sacramento.gpkg"), delete_dsn = TRUE, quiet = TRUE)

logmsg("CPAD holdings within study area + 2 km: ", nrow(cpad), " | distinct park units (UNIT_ID): ", length(unique(cpad$UNIT_ID)))
logmsg("Holdings with YR_EST filled: ", sum(!is.na(cpad$YR_EST) & cpad$YR_EST > 0), " of ", nrow(cpad))
cat("\nLAND_WATER:\n"); print(table(cpad$LAND_WATER, useNA = "ifany"))
cat("\nAgency level:\n"); print(table(cpad$AGNCY_LEV, useNA = "ifany"))
cat("\nAccess type:\n"); print(table(cpad$ACCESS_TYP, useNA = "ifany"))
cat("\nHolding area (ha) quantiles:\n"); print(round(quantile(cpad$area_ha, c(0, .1, .25, .5, .75, .9, 1)), 2))
logmsg("STAGE_01_DONE")
