#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/conus/05_terraclimate_deficit.R   (variant "conus", step 5 of 9)
##
## WHAT   Pulls the monthly TerraClimate climatic water deficit (def, mm) and
##        vapour-pressure deficit (vpd, kPa) 1958-latest at every cohort plot
##        through the OPeNDAP service (nearest 4 km cell; windowed reads over
##        the plots' bounding box, one per variable and year) and caches the
##        monthly series. Step 6 turns the deficit series into the yearly
##        growing-season anomaly cwd_z.
## WHY    Port of FINN-fia `scripts/03c_attach_drought.R --full` (the cache is
##        what the chain needs; the per-interval drought summaries the original
##        also wrote are not used by the annual preparation). TerraClimate is a
##        live API: the cache file is the reproducible layer, its sha256 is in
##        PROVENANCE.md. The cache is reused when it covers exactly the same
##        plots, years and variables; otherwise everything is pulled again
##        (~130 reads, 1-3 h depending on the server).
## INPUT  <DATA_DIR>/intermediate/conus/calibration_set_plots.csv (step 2);
##        TerraClimate via https://thredds.northwestknowledge.net (PROVENANCE.md)
## OUTPUT <DATA_DIR>/intermediate/conus/drought_monthly_cache.rds
##          list(plotID, years, vars, monthly = list(def = [plots x 12 x years], vpd = ...))
## RUN    Rscript calibration/conus/05_terraclimate_deficit.R   (from the repository root)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages({ library(data.table); library(ncdf4) })

DATA_DIR <- {
  cands <- c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))
  d <- cands[dir.exists(cands)][1]
  if (is.na(d)) stop("FINN-zoo-data/fia not found (expected two levels above the repository, see calibration/README.md)")
  normalizePath(d)
}
INT <- file.path(DATA_DIR, "intermediate", "conus")

PARAMS <- list(
  vars          = c(def = "def", vpd = "vpd"),   # TerraClimate variable names: climatic water deficit, vapour-pressure deficit
  baseline_from = 1958L,                         # TerraClimate begins 1958: the plot climatology for the z-anomaly
  baseline_to   = 2023L,                         # clamped to the latest year the server has
  base_url      = "https://thredds.northwestknowledge.net/thredds/dodsC/TERRACLIMATE_ALL/data",
  cache_rds     = file.path(INT, "drought_monthly_cache.rds")
)

plots <- fread(file.path(INT, "calibration_set_plots.csv"))
plots <- plots[is.finite(lat) & is.finite(lon) & is.finite(first_year) & is.finite(last_year)]
message(sprintf("[plots] %d plots; lon [%.2f, %.2f], lat [%.2f, %.2f]", nrow(plots), min(plots$lon), max(plots$lon), min(plots$lat), max(plots$lat)))

## reuse the cache when it covers exactly these plots (the cache written for the
## annual-v5 cohort is the reproducible layer of this live source)
if (file.exists(PARAMS$cache_rds)) {
  cache <- readRDS(PARAMS$cache_rds)
  if (identical(cache$plotID, plots$plotID) && identical(cache$vars, PARAMS$vars) &&
      min(cache$years) == PARAMS$baseline_from && max(cache$years) <= PARAMS$baseline_to) {
    message(sprintf("[cache] HIT: %s covers these %d plots, years %d..%d; nothing to pull",
                    PARAMS$cache_rds, length(cache$plotID), min(cache$years), max(cache$years)))
    quit(save = "no", status = 0)
  }
  message("[cache] present but for a different plot set or years; pulling again")
}

latest_year <- PARAMS$baseline_to
for (y in PARAMS$baseline_to:(PARAMS$baseline_to - 3L)) {
  u <- sprintf("%s/TerraClimate_%s_%d.nc", PARAMS$base_url, PARAMS$vars[1], y)
  ok <- tryCatch({ nc <- nc_open(u); nc_close(nc); TRUE }, error = function(e) FALSE)
  if (ok) { latest_year <- y; break }
}
years <- PARAMS$baseline_from:latest_year
message(sprintf("[terraclimate] latest available year = %d; pulling %d..%d", latest_year, min(years), max(years)))

bbox <- list(lon = range(plots$lon) + c(-0.1, 0.1), lat = range(plots$lat) + c(-0.1, 0.1))
GRID <- local({
  u <- sprintf("%s/TerraClimate_%s_%d.nc", PARAMS$base_url, PARAMS$vars[1], latest_year)
  nc <- nc_open(u); on.exit(nc_close(nc))
  lonv <- nc$dim$lon$vals; latv <- nc$dim$lat$vals
  li <- which(lonv >= bbox$lon[1] & lonv <= bbox$lon[2]); la <- which(latv >= bbox$lat[1] & latv <= bbox$lat[2])
  lon_sub <- lonv[li]; lat_sub <- latv[la]
  list(lon0 = min(li), lat0 = min(la), nlon = length(li), nlat = length(la),
       pil = vapply(plots$lon, function(x) which.min(abs(lon_sub - x)), integer(1)),
       pjl = vapply(plots$lat, function(y) which.min(abs(lat_sub - y)), integer(1)))
})
extract_var_year <- function(var, year) {
  u <- sprintf("%s/TerraClimate_%s_%d.nc", PARAMS$base_url, var, year)
  nc <- nc_open(u); on.exit(nc_close(nc))
  arr <- ncvar_get(nc, var, start = c(GRID$lon0, GRID$lat0, 1), count = c(GRID$nlon, GRID$nlat, 12))
  res <- matrix(NA_real_, nrow(plots), 12L)
  for (m in 1:12) res[, m] <- arr[cbind(GRID$pil, GRID$pjl, m)]
  res
}
monthly <- list()
for (v in names(PARAMS$vars)) {
  arr <- array(NA_real_, dim = c(nrow(plots), 12L, length(years)), dimnames = list(NULL, NULL, as.character(years)))
  for (iy in seq_along(years)) {
    arr[, , iy] <- tryCatch(extract_var_year(PARAMS$vars[v], years[iy]),
                            error = function(e) { message("  ! ", v, " ", years[iy], ": ", conditionMessage(e)); NA_real_ })
    if (iy %% 10 == 0) message(sprintf("  %s: %d/%d years", v, iy, length(years)))
  }
  monthly[[v]] <- arr
}
saveRDS(list(plotID = plots$plotID, years = years, vars = PARAMS$vars, monthly = monthly), PARAMS$cache_rds)
message("[out] ", PARAMS$cache_rds, " (pulled ", format(Sys.Date()), ")")
