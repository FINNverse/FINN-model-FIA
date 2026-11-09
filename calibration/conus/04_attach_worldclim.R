#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/conus/04_attach_worldclim.R   (variant "conus", step 4 of 9)
##
## WHAT   Attaches the WorldClim 2.1 normals (1970-2000, 30 arc-sec) to every
##        plot of the metadata table: tavg_c (BIO1, annual mean temperature),
##        prec_mm (BIO12, annual precipitation), and carries the plot
##        coordinates as x (= lon) and y (= lat).
## WHY    Port of FINN-fia `scripts/03_attach_environment.R`, block A only.
##        In the annual preparation (step 8) these static normals are used to
##        drop plots outside a WorldClim cell and to carry x, y; the yearly
##        climate the model sees comes from TerraClimate (steps 5-6). The
##        ecoregion block of the original is not needed by the FINN tables and
##        is not ported; elevation is not used either.
## INPUT  <DATA_DIR>/intermediate/conus/plot_metadata.csv (step 1);
##        WorldClim 2.1 rasters wc2.1_30s_bio_1.tif, wc2.1_30s_bio_12.tif
##        (PROVENANCE.md), under <DATA_DIR>/raw/worldclim or a `worldclim`
##        folder next to the models folder (see WORLDCLIM_DIR below).
## OUTPUT <DATA_DIR>/intermediate/conus/plot_env.csv   plotID, x, y, tavg_c, prec_mm
## RUN    Rscript calibration/conus/04_attach_worldclim.R   (from the repository root, ~30 s)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages({ library(data.table); library(terra) })

DATA_DIR <- {
  cands <- c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))
  d <- cands[dir.exists(cands)][1]
  if (is.na(d)) stop("FINN-zoo-data/fia not found (expected two levels above the repository, see calibration/README.md)")
  normalizePath(d)
}
WORLDCLIM_DIR <- {
  cands <- c(file.path(DATA_DIR, "raw", "worldclim"), file.path("..", "..", "worldclim"))
  d <- cands[dir.exists(cands)][1]
  if (is.na(d)) stop("WorldClim folder not found under ", paste(cands, collapse = " or "))
  normalizePath(d)
}
INT <- file.path(DATA_DIR, "intermediate", "conus")

PARAMS <- list(
  rasters = c(tavg_c = "bio_1", prec_mm = "bio_12"),   # WorldClim 2.1 bioclim: annual mean temperature (C), annual precipitation (mm)
  resolution = "30s"                                    # 30 arc-seconds (~1 km), the finest WorldClim product
)
## wc2.1_<res>_<stem>.tif, either directly in the folder or in its wc2.1_<res>_bio/ subfolder
find_raster <- function(stem) {
  direct <- file.path(WORLDCLIM_DIR, sprintf("wc2.1_%s_%s.tif", PARAMS$resolution, stem))
  subdir <- file.path(WORLDCLIM_DIR, sprintf("wc2.1_%s_bio", PARAMS$resolution), sprintf("wc2.1_%s_%s.tif", PARAMS$resolution, stem))
  p <- c(subdir, direct)[file.exists(c(subdir, direct))][1]
  if (is.na(p)) stop("WorldClim raster not found: ", stem, " at ", PARAMS$resolution)
  p
}

md <- fread(file.path(INT, "plot_metadata.csv"), select = c("plotID", "lat", "lon"))[!is.na(lat) & !is.na(lon)]
pts <- terra::vect(md, geom = c("lon", "lat"), crs = "EPSG:4326")
for (v in names(PARAMS$rasters)) {
  p <- find_raster(PARAMS$rasters[[v]])
  message("[read] ", v, " <- ", basename(p))
  md[, (v) := terra::extract(terra::rast(p), pts)[, 2]]
}
md[, `:=`(x = lon, y = lat)]
plot_env <- md[, .(plotID, x, y, tavg_c, prec_mm)]
fwrite(plot_env, file.path(INT, "plot_env.csv"))
message(sprintf("[out] plot_env.csv: %d plots; tavg_c [%.1f, %.1f] (NA %d), prec_mm [%.0f, %.0f] (NA %d)",
                nrow(plot_env), min(plot_env$tavg_c, na.rm = TRUE), max(plot_env$tavg_c, na.rm = TRUE), sum(is.na(plot_env$tavg_c)),
                min(plot_env$prec_mm, na.rm = TRUE), max(plot_env$prec_mm, na.rm = TRUE), sum(is.na(plot_env$prec_mm))))
