#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/conus/06_terraclimate_annual.R   (variant "conus", step 6 of 9)
##
## WHAT   Builds the yearly climate of every cohort plot and calendar year
##        2000-latest: tavg_c (mean over months of (tmax + tmin) / 2), prec_mm
##        (annual sum of ppt) from TerraClimate monthly tmax / tmin / ppt pulled
##        through OPeNDAP (cached), and cwd_z, the growing-season (Apr-Sep)
##        climatic water deficit of the year as a z-anomaly against the plot's
##        1958-latest climatology, from the deficit cache of step 5.
## WHY    Port of FINN-fia `scripts/03d_attach_annual_climate.R`. These are the
##        three covariates of the conus model.
## INPUT  <DATA_DIR>/intermediate/conus/calibration_set_plots.csv (step 2),
##        drought_monthly_cache.rds (step 5); TerraClimate via
##        https://thredds.northwestknowledge.net (PROVENANCE.md)
## OUTPUT <DATA_DIR>/intermediate/conus/annual_climate_cache.rds   monthly tmax/tmin/ppt [plots x 12 x years]
##        <DATA_DIR>/intermediate/conus/annual_climate_by_plot.csv  plotID, year, tavg_c, prec_mm, cwd_z, cwd_severe
## RUN    Rscript calibration/conus/06_terraclimate_annual.R   (from the repository root;
##        seconds with the cache, ~1 h without: 75 OPeNDAP reads)
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
  vars      = c(tmax = "tmax", tmin = "tmin", ppt = "ppt"),   # TerraClimate monthly maximum / minimum temperature, precipitation
  year_from = 2000L,            # first inventory year of the cohort
  year_to   = 2024L,            # clamped to the latest year the server has
  gs_months = 4:9,              # growing season for the deficit anomaly
  sev_z     = 1.0,              # cwd_severe flag threshold (carried along, not a model covariate)
  base_url  = "https://thredds.northwestknowledge.net/thredds/dodsC/TERRACLIMATE_ALL/data",
  cache_rds = file.path(INT, "annual_climate_cache.rds"),
  def_cache = file.path(INT, "drought_monthly_cache.rds"),
  out_csv   = file.path(INT, "annual_climate_by_plot.csv")
)

plots <- fread(file.path(INT, "calibration_set_plots.csv"))[is.finite(lat) & is.finite(lon)]
message(sprintf("[plots] %d plots", nrow(plots)))

## ---- monthly tmax / tmin / ppt: cache or OPeNDAP ---------------------------
cache_ok <- FALSE
if (file.exists(PARAMS$cache_rds)) {
  cache <- readRDS(PARAMS$cache_rds)
  cache_ok <- identical(cache$plotID, plots$plotID) && identical(cache$vars, PARAMS$vars) &&
              min(cache$years) == PARAMS$year_from && max(cache$years) <= PARAMS$year_to
}
if (cache_ok) {
  message(sprintf("[cache] HIT: %s (years %d..%d)", PARAMS$cache_rds, min(cache$years), max(cache$years)))
  monthly <- cache$monthly; years <- cache$years
} else {
  latest_year <- PARAMS$year_to
  for (y in PARAMS$year_to:(PARAMS$year_to - 3L)) {
    u <- sprintf("%s/TerraClimate_%s_%d.nc", PARAMS$base_url, PARAMS$vars[1], y)
    ok <- tryCatch({ nc <- nc_open(u); nc_close(nc); TRUE }, error = function(e) FALSE)
    if (ok) { latest_year <- y; break }
  }
  years <- PARAMS$year_from:latest_year
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
      if (iy %% 5 == 0) message(sprintf("  %s: %d/%d years", v, iy, length(years)))
    }
    monthly[[v]] <- arr
  }
  saveRDS(list(plotID = plots$plotID, years = years, vars = PARAMS$vars, monthly = monthly), PARAMS$cache_rds)
  message("[out] ", PARAMS$cache_rds, " (pulled ", format(Sys.Date()), ")")
}

## ---- annual temperature and precipitation ----------------------------------
tavg <- apply((monthly$tmax + monthly$tmin) / 2, c(1, 3), mean, na.rm = TRUE)
prec <- apply(monthly$ppt, c(1, 3), sum, na.rm = TRUE)
prec[apply(is.na(monthly$ppt), c(1, 3), all)] <- NA_real_

## ---- yearly deficit anomaly from the step-5 cache --------------------------
if (!file.exists(PARAMS$def_cache)) stop("deficit cache not found: ", PARAMS$def_cache, " -- run 05_terraclimate_deficit.R first")
dc <- readRDS(PARAMS$def_cache)
idx <- match(plots$plotID, dc$plotID)
if (anyNA(idx)) stop(sum(is.na(idx)), " plots are not in the deficit cache -- re-run 05_terraclimate_deficit.R on this cohort")
def_gs <- apply(dc$monthly$def[idx, PARAMS$gs_months, , drop = FALSE], c(1, 3), sum, na.rm = TRUE)
mu <- rowMeans(def_gs, na.rm = TRUE); sdv <- apply(def_gs, 1, sd, na.rm = TRUE)
sdv[!is.finite(sdv) | sdv == 0] <- NA_real_
zdef <- (def_gs - mu) / sdv
colnames(zdef) <- as.character(dc$years)
cwd_z <- matrix(NA_real_, nrow(plots), length(years), dimnames = list(NULL, as.character(years)))
have <- intersect(as.character(years), colnames(zdef))
cwd_z[, have] <- zdef[, have, drop = FALSE]

ann <- rbindlist(lapply(seq_along(years), function(iy) data.table(
  plotID = plots$plotID, year = years[iy], tavg_c = tavg[, iy], prec_mm = prec[, iy], cwd_z = cwd_z[, iy])))
ann[, cwd_severe := as.integer(cwd_z > PARAMS$sev_z)]
for (cc in c("tavg_c", "prec_mm", "cwd_z")) set(ann, which(!is.finite(ann[[cc]])), cc, NA_real_)
setorder(ann, plotID, year)
fwrite(ann, PARAMS$out_csv)
message(sprintf("[out] %s: %d plots x %d years; NA: tavg %d, prec %d, cwd_z %d; severe-deficit plot-years %.1f %%",
                PARAMS$out_csv, uniqueN(ann$plotID), length(years), sum(is.na(ann$tavg_c)), sum(is.na(ann$prec_mm)),
                sum(is.na(ann$cwd_z)), 100 * mean(ann$cwd_severe, na.rm = TRUE)))
