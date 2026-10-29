# _common.R  -- shared settings sourced by every stage script
lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", "4.6")
.libPaths(c(lib, .libPaths()))
ROOT     <- Sys.getenv("PCA_ROOT", unset = "path/to/your/project/folder")   # where data_raw/, data_proc/, outputs/ and logs/ live
DIR_RAW  <- file.path(ROOT, "data_raw")
DIR_PROC <- file.path(ROOT, "data_proc")
DIR_OUT  <- file.path(ROOT, "outputs")
DIR_LOG  <- file.path(ROOT, "logs")
SEED     <- 20261006L                      # recorded seed for the whole pipeline
CRS_WORK <- "EPSG:32610"                   # WGS 84 / UTM 10N = native Landsat grid for Sacramento
CBD_LONLAT <- c(-121.49333, 38.57667)      # CBD defined as the California State Capitol
options(timeout = 3600, warn = 1)
suppressPackageStartupMessages({ library(sf); library(terra); library(data.table); library(jsonlite) })
logmsg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), paste0(...), "\n")
stopif <- function(cond, msg) if (isTRUE(cond)) stop("VALIDATION FAILED: ", msg, call. = FALSE)
