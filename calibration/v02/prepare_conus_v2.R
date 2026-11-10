#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/prepare_conus_v2.R   (variant "conus", data v2, step B)
##
## V0.2 DATA (2026-11-10): copy of calibration/conus/09_prepare_conus.R reading
##        `annual-v6-ingrowth` (calibration/v02/prepare_conus_v2_annual.R, the
##        complete ingrowth target) instead of `annual-v5`, and writing
##        `conus-annual-v2`. Nothing else changes. The step-9 header follows.
##
## WHAT   Turns the prepared set `annual-v5` written by step 8 (conterminous
##        US, 20 species, annual time step, per-site remeasurement intervals,
##        spatially blocked train / holdout) into the prepared set of this
##        variant: `conus-annual-v1`. Nothing is re-selected here; the script
##        re-shapes, drops what the model does not use, checks the invariants
##        FINN relies on and writes the summary a rerun is diffed against.
## WHY    annual-v5 is the data of the FINN-fia FIA fits (2026-09); steps 1-8
##        of this folder rebuild it from the public sources (PROVENANCE.md).
##        What changes for the tool:
##        - the environment is written in NATURAL units (the FINN-fia fits used
##          the pre-standardised `env_dt.csv`); FINN standardises internally
##          (`env_autoscale`) and stores the constants with the model, so users
##          pass natural units;
##        - the plot coordinates `x`, `y` are NOT model covariates (the FINN-fia
##          fits had them as network inputs and predicted the holdout with them
##          held at the training mean; for a tool a pure climate model is
##          simpler to apply and to transfer). They stay in siteID_dt.csv for maps.
## INPUT  <DATA_DIR>/prepared/annual-v5/{train1500,holdout_blocked}/ (step 8), plus
##        site_blocks.csv, holdout_separation.csv and PREPARATION.md in <DATA_DIR>/prepared/annual-v5/.
## OUTPUT <DATA_DIR>/prepared/conus-annual-v1/{train,holdout}/
##          obs_dt.csv          site x census x species (year = simulated year of the census,
##                              period_length = years since the previous census): ba, dbh,
##                              trees, growth (annual relative dbh increment), growth_n,
##                              mort (annual rate), n_at_risk, n_died, reg (recruits per ha
##                              per interval), species_name
##          tree_obs_dt.csv     one row per tree alive at the start of an interval: siteID, year
##                              (census closing the interval), patchID, species, dbh (start),
##                              growth (annual relative), died (0/1), period_length, trees
##          env_dt.csv          site x year (1..last census), natural units: tavg_c, prec_mm, cwd_z
##          init_cohort_dt.csv  trees of the first inventory (siteID, patchID, species, dbh, trees, treeName)
##          species_dt.csv, siteID_dt.csv (siteID <-> FIA plot id, lon, lat), site_calendar.csv
##        <DATA_DIR>/prepared/conus-annual-v1/PREPARATION.md, site_blocks.csv, holdout_separation.csv
##        calibration/conus/expected_preparation.csv
## RUN    Rscript calibration/conus/09_prepare_conus.R      (from the repository root, ~20 s)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({
  library(FINN)
  library(data.table)
})

DATA_DIR <- {
  cands <- c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))
  d <- cands[dir.exists(cands)][1]
  if (is.na(d)) stop("FINN-zoo-data/fia not found (expected two levels above the repository, see calibration/README.md)")
  normalizePath(d)
}

PARAMS <- list(
  source_version = "annual-v6-ingrowth",    # v0.2: step 8 with the complete ingrowth target
  prep_version   = "conus-annual-v2",
  sets           = c(train = "train1500", holdout = "holdout_blocked"),
  ## the model's covariates (yearly TerraClimate per plot): mean annual
  ## temperature, annual precipitation, growing-season water-deficit z-anomaly
  env_vars       = c("tavg_c", "prec_mm", "cwd_z"),
  patch_size_ha  = 0.01681,                 # one FIA subplot (24 ft radius)
  n_patches      = 4L
)
SRC <- file.path(DATA_DIR, "prepared", PARAMS$source_version)
OUT <- file.path(DATA_DIR, "prepared", PARAMS$prep_version)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

species <- fread(file.path(SRC, PARAMS$sets[["train"]], "species_dt.csv"))[order(species)]
stopifnot(identical(species$species, seq_len(nrow(species))))

summaries <- list()
for (nm in names(PARAMS$sets)) {
  s <- file.path(SRC, PARAMS$sets[[nm]]); d <- file.path(OUT, nm); dir.create(d, showWarnings = FALSE)
  obs   <- fread(file.path(s, "obs_dt.csv"))
  envu  <- fread(file.path(s, "env_unscaled_dt.csv"))
  trees <- fread(file.path(s, "tree_obs_dt.csv"))
  init  <- fread(file.path(s, "init_cohort_dt.csv"))
  cal   <- fread(file.path(s, "site_calendar.csv"))
  ids   <- fread(file.path(s, "siteID_dt.csv"))
  stopifnot(identical(fread(file.path(s, "species_dt.csv"))[order(species)], species))

  ## invariants FINN relies on
  stopifnot(all(obs$species %in% species$species),
            all(obs[, paste(siteID, year)] %in% envu[, paste(siteID, year)]),   # env covers every observation
            all(trees[, paste(siteID, year)] %in% obs[, paste(siteID, year)]),  # every tree interval has its census row
            all(init$siteID %in% obs$siteID), all(trees$died %in% c(0, 1)),
            all(obs$n_died <= obs$n_at_risk), !anyNA(obs$period_length))
  per_site_env <- envu[, .(n_years = .N, first = min(year), last = max(year)), by = siteID]
  stopifnot(all(per_site_env$first == 1L), all(per_site_env$n_years == per_site_env$last))

  o <- obs[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, species_name, period_length)]
  setorder(o, siteID, year, species)
  e <- envu[, c("siteID", "year", PARAMS$env_vars), with = FALSE]; setorder(e, siteID, year)
  tr <- trees[, .(siteID, year, patchID, species, dbh, growth, died, period_length, trees)]
  i  <- init[, .(siteID, patchID, species, dbh, trees, treeName, species_name, living, status)]
  site_key <- unique(ids[, .(siteID, siteName)])
  site_key <- merge(site_key, unique(envu[, .(siteID, lon = x, lat = y)]), by = "siteID")
  setorder(site_key, siteID)

  fwrite(o, file.path(d, "obs_dt.csv")); fwrite(e, file.path(d, "env_dt.csv"))
  fwrite(tr, file.path(d, "tree_obs_dt.csv")); fwrite(i, file.path(d, "init_cohort_dt.csv"))
  fwrite(species, file.path(d, "species_dt.csv")); fwrite(site_key, file.path(d, "siteID_dt.csv"))
  fwrite(cal, file.path(d, "site_calendar.csv"))
  message(sprintf("%-8s %d sites, %d census rows (%d site x census intervals, period lengths %s), %d tree intervals, %d init trees, env %d..%d years per site",
                  nm, uniqueN(o$siteID), nrow(o), uniqueN(o[, paste(siteID, year)]),
                  paste(names(table(unique(o[, .(siteID, year, period_length)])$period_length)), collapse = "/"),
                  nrow(tr), nrow(i), min(per_site_env$n_years), max(per_site_env$n_years)))

  ## the summary a rerun elsewhere is diffed against
  o <- merge(o, cal[, .(siteID, year = rel_year, period)], by = c("siteID", "year"))
  stopifnot(nrow(o) == nrow(obs))
  tot <- o[, .(set = nm, level = "total", species = "all", period = "all",
               n_sites = uniqueN(siteID), n_rows = .N,
               sum_ba = round(sum(ba, na.rm = TRUE), 4), sum_trees = round(sum(trees, na.rm = TRUE), 1),
               sum_reg = round(sum(reg, na.rm = TRUE), 2),
               n_at_risk = sum(n_at_risk), n_died = sum(n_died), growth_n = sum(growth_n),
               mean_growth = round(weighted.mean(growth, growth_n, na.rm = TRUE), 6))]
  sp <- o[, .(level = "species", n_sites = uniqueN(siteID[trees > 0]), n_rows = .N,
              sum_ba = round(sum(ba, na.rm = TRUE), 4), sum_trees = round(sum(trees, na.rm = TRUE), 1),
              sum_reg = round(sum(reg, na.rm = TRUE), 2),
              n_at_risk = sum(n_at_risk), n_died = sum(n_died), growth_n = sum(growth_n),
              mean_growth = round(weighted.mean(growth, growth_n, na.rm = TRUE), 6)),
          by = .(species = species_name, period = as.character(period - 1L))][, set := nm]
  trs <- tr[, .(set = nm, level = "tree_intervals", species = "all", period = "all", n_sites = uniqueN(siteID), n_rows = .N,
                sum_ba = round(sum(pi * (dbh / 200)^2), 4), sum_trees = sum(trees), sum_reg = NA_real_,
                n_at_risk = .N, n_died = sum(died), growth_n = sum(is.finite(growth)),
                mean_growth = round(mean(growth, na.rm = TRUE), 6))]
  ini <- i[living == TRUE, .(set = nm, level = "init_trees", species = "all", period = "0", n_sites = uniqueN(siteID), n_rows = .N,
                             sum_ba = round(sum(pi * (dbh / 200)^2), 4), sum_trees = sum(trees), sum_reg = NA_real_,
                             n_at_risk = NA_integer_, n_died = NA_integer_, growth_n = NA_integer_, mean_growth = NA_real_)]
  summaries[[nm]] <- rbind(tot, sp[, names(tot), with = FALSE], trs, ini)
}
expected <- rbindlist(summaries)[, variant := "conus_v2"]
fwrite(expected, file.path("calibration", "v02", "expected_preparation_conus_v2.csv"))

## carried over unchanged: blocks and separation of the two sets
file.copy(file.path(SRC, c("site_blocks.csv", "holdout_separation.csv")), OUT, overwrite = TRUE)
sep <- fread(file.path(SRC, "holdout_separation.csv"))

writeLines(c(
  sprintf("# Prepared data: %s", PARAMS$prep_version), "",
  sprintf("Written %s by `calibration/conus/09_prepare_conus.R` (FINN %s) from the prepared set `%s` of step 8.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), as.character(packageVersion("FINN")), PARAMS$source_version), "",
  "## Source", "",
  "- `prepared/annual-v5` (train1500, holdout_blocked), written by `calibration/conus/08_prepare_annual_v5.R`, the port of",
  "  FINN-fia `scripts/07_prepare_finn_inputs.R` (commit 58fcf3e, 2026-09-10), from the FIADB CSV dump (USDA Forest Service FIA",
  "  DataMart, public domain / CC0, dump of June 2025), WorldClim 2.1 and yearly TerraClimate (1958-2023) summaries per plot",
  "  (calibration/conus steps 1-8, calibration/PROVENANCE.md). The PREPARATION.md of annual-v5 is appended below.",
  "- Cohort: FIA plots with >= 2 inventories, all four subplots consistently sampled, condition-level management tier",
  "  M0/M1, the 20 species of the FINN-fia species ranking jointly >= 99 % of the plot's basal area at the first",
  "  inventory (plots dominated by other species are excluded; the 20 species are therefore modelled without an",
  "  `other` class). Pacific-Northwest macroplot records thinned 1 in 6.", "",
  "## Definitions", "",
  "- Site = FIA plot, patch = one of its four subplots (0.01681 ha). Annual time step; simulated year 0 is the first",
  "  inventory, a census at simulated year k closes an interval of `period_length` years (5-10, mostly 10 in the West,",
  "  5-7 in the East); up to 25 years per site.",
  "- `growth`: mean annual relative dbh increment, (1 + dbh/dbh_before - 1)^(1/p) - 1; `mort`: annual rate",
  "  1 - (1 - fraction died)^(1/p) with the closed-cohort counts `n_at_risk`/`n_died` over the interval; `reg`: ingrowth",
  "  into the 12.7 cm tally threshold per ha per interval. Tree table: dbh at the interval start, annual relative growth,",
  "  died 0/1 over the interval.",
  sprintf("- Environment, natural units, yearly per plot: `%s`. tavg_c = TerraClimate mean annual temperature of the year (C),",
          paste(PARAMS$env_vars, collapse = "`, `")),
  "  prec_mm = annual precipitation (mm), cwd_z = z-anomaly of the growing-season climatic water deficit of the year relative",
  "  to the plot's 1958-2023 distribution (> 1 = severe drought year). Plot coordinates are NOT covariates of this variant.", "",
  "## Sets", "",
  sprintf("- train: %d sites (FINN-fia train1500, seed 1, drawn from the retained 2-degree blocks)", expected[set == "train" & level == "total", n_sites]),
  sprintf("- holdout: %d sites (FINN-fia holdout_blocked, seed 11, drawn from the 26 of 128 held-out 2-degree blocks, block seed 21)",
          expected[set == "holdout" & level == "total", n_sites]), "",
  sprintf("Distance from each holdout site to the nearest training site: median %.0f km, minimum %.0f km, 10 %% below %.0f km.",
          median(sep$nearest_train_km), min(sep$nearest_train_km), quantile(sep$nearest_train_km, 0.1)), "",
  "## PARAMS", "", "```r", capture.output(str(PARAMS)), "```", "",
  "---", "", "# Appendix: PREPARATION.md of annual-v5 (step 8)", "",
  readLines(file.path(SRC, "PREPARATION.md"))
), file.path(OUT, "PREPARATION.md"))
message("done: ", OUT)
