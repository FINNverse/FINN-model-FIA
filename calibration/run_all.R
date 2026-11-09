#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/run_all.R
##
## WHAT   Runs the whole calibration chain of this release in its fixed order,
##        every step as its own `Rscript`, so each script stays self-contained
##        and this file is nothing but the order:
##          1. data preparation of the oregon variant (1 step) and of the conus
##             variant (9 steps: FIADB metadata -> cohort -> trees -> WorldClim ->
##             TerraClimate deficit -> yearly climate -> species ranking ->
##             annual-v5 -> conus-annual-v1)
##          2. the production fit of each variant (one fit per variant, no repeats)
##          3. the validation report of each variant
##          4. the package files (inst/model) and the combined expected_preparation.csv
##        A step that fails stops the chain. Scripts that are already done can
##        be skipped by setting `skip` below; nothing else is configurable here.
## NEEDS  the raw layers listed in calibration/PROVENANCE.md (FIADB dump,
##        WorldClim 2.1 rasters, the Oregon export; TerraClimate is pulled or
##        read from its cache) under ../../FINN-zoo-data/fia relative to the
##        repository root, and FINN 0.2.0 (commit 0848b60) in ~/R/lib-finn-0848b60.
## RUN    Rscript calibration/run_all.R         (from the repository root)
##        Timing, single-threaded: conus preparation ~30 min with the
##        TerraClimate caches (hours without, ~200 OPeNDAP reads), oregon
##        preparation ~1 min, oregon fit ~10 min, conus fit ~3.5 h, oregon
##        validation ~5 min, conus validation ~30 min.
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")   # technical: every step is single-threaded by design
steps <- c(
  "calibration/oregon/01_prepare.R",
  "calibration/conus/01_fiadb_metadata.R",
  "calibration/conus/02_select_cohort.R",
  "calibration/conus/03_extract_trees.R",
  "calibration/conus/04_attach_worldclim.R",
  "calibration/conus/05_terraclimate_deficit.R",
  "calibration/conus/06_terraclimate_annual.R",
  "calibration/conus/07_species_ranking.R",
  "calibration/conus/08_prepare_annual_v5.R",
  "calibration/conus/09_prepare_conus.R",
  "calibration/oregon/02_fit.R",
  "calibration/conus/02_fit.R",
  "calibration/oregon/03_validation.R",
  "calibration/conus/03_validation.R",
  "calibration/04_package_files.R"
)
skip <- character(0)      # e.g. c("calibration/conus/02_fit.R") when the fit ran on a cluster

if (!file.exists("DESCRIPTION") || !dir.exists("calibration"))
  stop("run from the repository root of FINN-model-FIA")
for (s in steps) {
  if (s %in% skip) { message("skip  ", s); next }
  message(sprintf("[%s] start %s", format(Sys.time(), "%H:%M:%S"), s))
  status <- system2("Rscript", shQuote(s))
  if (status != 0) stop("step failed: ", s, " (exit status ", status, ")")
  message(sprintf("[%s] done  %s", format(Sys.time(), "%H:%M:%S"), s))
}
message("calibration chain complete")
