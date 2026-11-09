#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/04_package_files.R
##
## WHAT   Writes the files the installed package ships from the outputs of the
##        calibration chain, so the package and the calibration cannot drift:
##          inst/model/species.csv      variant, species index, name
##          inst/model/covariates.csv   variant, covariate, unit, source, derivation,
##                                      and the standardisation (center, scale) the
##                                      fit learned (env_autoscale), read by build_model()
##          inst/model/weights.json     per variant: release tag, parameter file, url,
##                                      sha256, FINN commit, patches, patch_size
##          inst/extdata/<variant>_example_{env,init}.csv   a few sites for the quickstart
##          calibration/expected_preparation.csv   both variants' preparation summaries
##        The release tag and the download url are placeholders until the
##        coordinator publishes the release; everything else is final.
## INPUT  calibration/<variant>/{env_scaling.csv, parameters_sha256.txt, expected_preparation.csv}
##        <DATA_DIR>/prepared/<version>/<train set>/species_dt.csv and the holdout set
##        (the example sites are holdout sites, which no fit has seen)
## RUN    Rscript calibration/04_package_files.R      (from the repository root, seconds)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

DATA_DIR <- {
  cands <- c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))
  d <- cands[dir.exists(cands)][1]
  if (is.na(d)) stop("FINN-zoo-data/fia not found (expected two levels above the repository, see calibration/README.md)")
  normalizePath(d)
}

PARAMS <- list(
  release_tag  = "v0.1.0",                 # placeholder until the coordinator tags the release
  url_template = "https://github.com/FINNverse/FINN-model-FIA/releases/download/%s/%s",
  finn_commit  = "0848b60",                # FINN 0.2.0 main, the version the parameters were fitted under
  threads      = 1L,                       # the fits and the quickstart run single-threaded (torch_set_num_threads)
  sim_batchsize = 250L,                    # sites per simulateForest() batch pinned in the quickstart (memory; draws depend on it)
  n_example_sites = 6L,                    # holdout sites shipped as the quickstart example
  variants = list(
    oregon = list(prep_version = "oregon-10yr-v1", train = "train_blocked", holdout = "holdout_blocked",
                  patches = 4L, patch_size = 0.01681, step_years = 10L,
                  covariates = data.table(
                    name = c("temp", "prec"), unit = c("degree C", "mm per year"),
                    source = c("WorldClim 2.1 bio1 (1970-2000 normal), 30 arc-sec", "WorldClim 2.1 bio12 (1970-2000 normal), 30 arc-sec"),
                    derivation = c("mean annual temperature at the plot; constant over the simulation steps",
                                   "annual precipitation sum at the plot; constant over the simulation steps"))),
    conus  = list(prep_version = "conus-annual-v1", train = "train", holdout = "holdout",
                  patches = 4L, patch_size = 0.01681, step_years = 1L,
                  covariates = data.table(
                    name = c("tavg_c", "prec_mm", "cwd_z"), unit = c("degree C", "mm per year", "z-score"),
                    source = c("TerraClimate (4 km), yearly", "TerraClimate (4 km), yearly", "TerraClimate (4 km), yearly"),
                    derivation = c("mean annual temperature of the simulated year at the plot",
                                   "annual precipitation sum of the simulated year at the plot",
                                   "growing-season (Apr-Sep) climatic water deficit of the year, standardised by the plot's 1958-2023 mean and sd; 0 = average year, > 1 = severe drought")))
  )
)
dir.create(file.path("inst", "model"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("inst", "extdata"), recursive = TRUE, showWarnings = FALSE)

species_all <- list(); cov_all <- list(); weights <- list(); expected <- list()
for (v in names(PARAMS$variants)) {
  p <- PARAMS$variants[[v]]
  train_dir <- file.path(DATA_DIR, "prepared", p$prep_version, p$train)
  hold_dir  <- file.path(DATA_DIR, "prepared", p$prep_version, p$holdout)
  sp <- fread(file.path(train_dir, "species_dt.csv"))[order(species)]
  species_all[[v]] <- data.table(variant = v, species = sp$species, name = sp$species_name)
  sc <- fread(file.path("calibration", v, "env_scaling.csv"))
  cov <- merge(p$covariates, sc[, .(name = variable, center, scale)], by = "name", sort = FALSE)
  stopifnot(nrow(cov) == nrow(p$covariates))
  cov_all[[v]] <- cbind(variant = v, cov)
  ## two release assets per variant: the name-matched parameter list (what
  ## load_model() rebuilds from) and the whole model object (torch_load() under
  ## the same FINN version restores it directly, methods included)
  sha <- readLines(file.path("calibration", v, "parameters_sha256.txt"))
  model_sha <- readLines(file.path("calibration", v, "model_sha256.txt"))
  file <- sprintf("parameters_%s.pt", v); model_file <- sprintf("model_%s.pt", v)
  weights[[v]] <- list(tag = PARAMS$release_tag, file = file, url = sprintf(PARAMS$url_template, PARAMS$release_tag, file),
                       sha256 = sha, model_file = model_file, model_url = sprintf(PARAMS$url_template, PARAMS$release_tag, model_file),
                       model_sha256 = model_sha, finn_commit = PARAMS$finn_commit, patches = p$patches, patch_size = p$patch_size,
                       step_years = p$step_years, data_version = p$prep_version,
                       threads = PARAMS$threads, sim_batchsize = PARAMS$sim_batchsize)
  file.copy(file.path(DATA_DIR, "fits", v, "parameters.pt"), file.path(DATA_DIR, "fits", v, file), overwrite = TRUE)
  file.copy(file.path(DATA_DIR, "fits", v, "model.pt"), file.path(DATA_DIR, "fits", v, model_file), overwrite = TRUE)
  expected[[v]] <- fread(file.path("calibration", v, "expected_preparation.csv"))
  ## example sites: the first holdout sites, env and initial trees, renumbered 1..n
  env <- fread(file.path(hold_dir, "env_dt.csv")); init <- fread(file.path(hold_dir, "init_cohort_dt.csv"))
  ids <- head(sort(unique(init$siteID)), PARAMS$n_example_sites)
  key <- data.table(siteID = ids, new = seq_along(ids))
  ex_env <- merge(env[siteID %in% ids], key, by = "siteID")[, siteID := new][, new := NULL]
  ex_init <- merge(init[siteID %in% ids], key, by = "siteID")[, siteID := new][, new := NULL]
  if ("living" %in% names(ex_init)) ex_init <- ex_init[living == TRUE]
  ex_init <- ex_init[, .(siteID, patchID, species, dbh, trees = if ("trees" %in% names(ex_init)) trees else 1L)]
  setorder(ex_env, siteID, year); setorder(ex_init, siteID, patchID, species)
  fwrite(ex_env, file.path("inst", "extdata", sprintf("%s_example_env.csv", v)))
  fwrite(ex_init, file.path("inst", "extdata", sprintf("%s_example_init.csv", v)))
}
fwrite(rbindlist(species_all), file.path("inst", "model", "species.csv"))
fwrite(rbindlist(cov_all), file.path("inst", "model", "covariates.csv"))
write_json(list(variants = weights), file.path("inst", "model", "weights.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
fwrite(rbindlist(expected, fill = TRUE), file.path("calibration", "expected_preparation.csv"))
message("inst/model and calibration/expected_preparation.csv written for variants: ", paste(names(PARAMS$variants), collapse = ", "))
