#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/run_all.R
##
## WHAT   Runs the whole calibration chain of this release in its fixed order,
##        every step as its own `Rscript`, so each script stays self-contained
##        and this file is nothing but the order:
##          1. data preparation of both variants (oregon, conus)
##          2. the production fit of each variant (one fit per variant, no repeats)
##          3. the validation report of each variant
##          4. the package files (inst/model) and the combined expected_preparation.csv
##        A step that fails stops the chain. Scripts that are already done can
##        be skipped by setting `skip` below; nothing else is configurable here.
## NEEDS  the raw layers listed in calibration/PROVENANCE.md under
##        ../../FINN-zoo-data/fia (relative to the repository root) and FINN
##        0.2.0 (commit 0848b60) in ~/R/lib-finn-0848b60.
## RUN    Rscript calibration/run_all.R         (from the repository root)
##        Timing, single-threaded: preparation ~2 min, oregon fit ~10 min,
##        conus fit ~3.5 h, oregon validation ~5 min, conus validation ~30 min.
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
steps <- c(
  "calibration/oregon/01_prepare.R",
  "calibration/conus/01_prepare.R",
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
