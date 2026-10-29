# 04b_park_classes.R
# Splits the screened (kept) CPAD units into the two treatment definitions and
# builds the control-exclusion layer. Runnable on its own.
#   Main:      urban parks  = >50% of unit area owned by City / County agency,
#              Recreation/Parks District, Community Services District or HOA,
#              and NOT golf course / cemetery / planned park.
#   Secondary: all kept green CPAD units (area >= 0.2 ha, NDVI > 0.3).
#   Control exclusion: ALL CPAD units (any status) + OSM green parks not in CPAD.
source("code/_common.R")   # run from the repository folder
set.seed(SEED)
h <- st_read(file.path(DIR_RAW, "cpad_holdings_sacramento.gpkg"), quiet = TRUE)
units <- st_read(file.path(DIR_PROC, "cpad_units_screened.gpkg"), quiet = TRUE)
ok_typ <- c("City Agency", "County Agency", "Recreation/Parks District", "Community Services District", "Home Owners Association")
bad_use <- "Golf Course|Cemetery|Planned Park"
hd <- as.data.table(st_drop_geometry(h))[, .(UNIT_ID, AGNCY_TYP, SPEC_USE, ha = area_ha)]
hd[, ok := AGNCY_TYP %in% ok_typ & !grepl(bad_use, SPEC_USE)]
share <- hd[, .(urban_share = sum(ha[ok]) / sum(ha), main_type = AGNCY_TYP[which.max(ha)],
                golf = any(grepl("Golf Course", SPEC_USE))), by = UNIT_ID]
units <- merge(units, share, by = "UNIT_ID", all.x = TRUE)
units$park_class <- ifelse(units$status != "kept", "not_treatment",
                     ifelse(units$urban_share > 0.5, "urban_park", "other_green_open_space"))
st_write(units, file.path(DIR_PROC, "cpad_units_classified.gpkg"), delete_dsn = TRUE, quiet = TRUE)
urban <- units[units$park_class == "urban_park", ]
allg  <- units[units$status == "kept", ]
st_write(urban, file.path(DIR_PROC, "parks_main_urban.gpkg"),  delete_dsn = TRUE, quiet = TRUE)
st_write(allg,  file.path(DIR_PROC, "parks_secondary_all_green.gpkg"), delete_dsn = TRUE, quiet = TRUE)

osm_only <- st_read(file.path(DIR_PROC, "osm_parks_not_in_cpad.gpkg"), quiet = TRUE)
excl <- c(st_geometry(units), st_geometry(osm_only))
excl <- st_sf(src = c(rep("CPAD", nrow(units)), rep("OSM_only", nrow(osm_only))), geometry = excl)
st_write(excl, file.path(DIR_PROC, "park_exclusion_all.gpkg"), delete_dsn = TRUE, quiet = TRUE)

cat("\nKept units by class:\n")
print(as.data.table(st_drop_geometry(allg))[, .(units = .N, area_ha = round(sum(area_ha)), median_ha = round(median(area_ha), 2),
       dated_after_1975 = sum(yr_est > 1975, na.rm = TRUE)), by = park_class])
cat("\nOther green open space, largest 10:\n")
og <- as.data.table(st_drop_geometry(allg))[park_class == "other_green_open_space"][order(-area_ha)][1:10, .(unit_name, main_type, golf, area_ha = round(area_ha))]
print(og)
cat("\nGolf-course units removed from main spec:", sum(units$golf & units$status == "kept", na.rm = TRUE), "\n")
stopif(nrow(urban) < 50, "fewer than 50 urban parks")
stopif(st_crs(urban)$epsg != 32610, "CRS lost")
logmsg("Urban parks: ", nrow(urban), " | all green: ", nrow(allg), " | exclusion polygons: ", nrow(excl))
logmsg("STAGE_04b_DONE")
