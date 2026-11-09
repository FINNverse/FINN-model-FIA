#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/oregon/02_fit.R   (variant "oregon")
##
## WHAT   Fits the Oregon FIA model: a fully mechanistic FINN (process functions
##        FINN::growth, FINN::mortality, FINN::regeneration_saturation, FINN::
##        competition) whose species parameters and linear climate responses
##        (mean annual temperature, annual precipitation) are learned from the
##        601 training plots of `oregon-10yr-v1`. One simulated step is a
##        10-year remeasurement interval.
## WHY    This is the FINN vignette set-up ("Fitting FINN to forest inventory
##        data", 200 sites) on the whole Oregon export, with three changes that
##        make it a tool rather than a demonstration: (1) the physical subplot
##        area as patch size after macroplot thinning (01_prepare.R), (2)
##        recruitment with a Beverton-Holt cap (regeneration_saturation) so a
##        bare-ground simulation cannot run away, (3) a spatially blocked
##        holdout. Climate enters linearly (`hidden = NULL`): with two
##        inventories per plot and 27 species, a per-species linear niche is
##        what the data can carry; the national variant has the networks.
## DATA   <DATA_DIR>/prepared/oregon-10yr-v1/train_blocked   (training)
##        <DATA_DIR>/prepared/oregon-10yr-v1/holdout_blocked (validation only, 03_validation.R)
## RUN    Rscript calibration/oregon/02_fit.R        (from the repository root; ~15 min single-threaded)
##        The cluster wrapper calibration/hpc/run_fit.sbatch passes only this path.
## OUTPUT <DATA_DIR>/fits/oregon/parameters.pt   named list of parameter tensors (the release asset)
##        <DATA_DIR>/fits/oregon/model.pt        the whole fitted object (validation convenience, not distributed)
##        <DATA_DIR>/fits/oregon/param_history.rds, loss_history.csv
##        calibration/oregon/fit_summary.md, loss_history.csv, env_scaling.csv, parameters_sha256.txt
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({
  library(FINN)
  library(data.table)
  library(torch)
})

DATA_DIR <- {
  cands <- c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))
  d <- cands[dir.exists(cands)][1]
  if (is.na(d)) stop("FINN-zoo-data/fia not found (expected two levels above the repository, see calibration/README.md)")
  normalizePath(d)
}
VARIANT <- "oregon"
FIT_DIR <- file.path(DATA_DIR, "fits", VARIANT); dir.create(FIT_DIR, recursive = TRUE, showWarnings = FALSE)
CAL_DIR <- file.path("calibration", VARIANT)

PARAMS <- list(
  prep_version  = "oregon-10yr-v1",
  train_set     = "train_blocked",
  holdout_set   = "holdout_blocked",       # not read here; recorded for the validation script
  patch_size_ha = 0.01681,                 # one FIA subplot, 24 ft radius
  n_patches     = 4L,                      # the four subplots of a plot
  recruits_dbh  = 12.9,                    # a recruit enters just above the 12.7 cm tally threshold
  env_formula   = ~ temp + prec,           # WorldClim mean annual temperature (C) and annual precipitation (mm)
  hidden        = NULL,                    # linear climate response per species (see WHY)
  ## Beverton-Holt recruitment cap, stems per ha per decade, initial value: the
  ## observed ingrowth is ~130 per ha per decade on average, so the cap starts
  ## well above it and is learned per species
  reg_K_init    = 500,
  ## recruitment a species gets without light (immigration), per ha per decade;
  ## FINN's default 0.2 is per year on annual data
  reg_floor     = 0.01,
  epochs        = 500L,                    # the vignette's length: every response plateaus by ~400 on this data
  batchsize     = 50L,                     # ~12 batches per epoch on 601 sites
  lr            = 0.01,
  ## every term is divided by its intercept-only baseline before the terms are
  ## summed, so the six responses start at equal weight (FINN >= 0.2.0 default);
  ## a numeric vector would multiply these baseline-scaled terms
  weights       = "auto",
  weight_scale  = "auto",
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse",
                    mortality = "binomial", regeneration = "nbinom"),   # FINN defaults, written out
  env_autoscale = TRUE,                    # natural units in, z-scaling learned and stored with the model
  device        = "cpu",
  threads       = 1L,                      # single-threaded: the fits share a node
  seed          = 42L
)
torch_set_num_threads(PARAMS$threads)

## ---- 1. data ---------------------------------------------------------------
train_dir <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs     <- fread(file.path(train_dir, "obs_dt.csv"))
env     <- fread(file.path(train_dir, "env_dt.csv"))
init    <- fread(file.path(train_dir, "init_cohort_dt.csv"))
species <- fread(file.path(train_dir, "species_dt.csv"))[order(species)]
n_species <- nrow(species)
init_cohort <- makeInitCohorts(init, Nspecies = n_species)
message(sprintf("train: %d sites, %d species, census steps %s, %d observation rows, %d initial trees",
                uniqueN(obs$siteID), n_species, paste(sort(unique(obs$year)), collapse = "/"), nrow(obs), nrow(init)))

## ---- 2. the model, written out ---------------------------------------------
FINN.seed(PARAMS$seed)
model <- finn(
  N_species    = n_species,
  recruits_dbh = PARAMS$recruits_dbh,
  reg_floor    = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(PARAMS$env_formula, FINN::growth,
                                       hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(PARAMS$env_formula, FINN::mortality,
                                       hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(PARAMS$env_formula, FINN::regeneration_saturation,
                                       hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       custom_parameters = list(reg_logK = rep(log(PARAMS$reg_K_init), n_species)))
)

## ---- 3. fit ----------------------------------------------------------------
t0 <- Sys.time()
fit(model,
    env           = env,
    data          = obs,
    init_cohort   = init_cohort,
    epochs        = PARAMS$epochs,
    batchsize     = PARAMS$batchsize,
    lr            = PARAMS$lr,
    weights       = PARAMS$weights,
    weight_scale  = PARAMS$weight_scale,
    loss_family   = PARAMS$loss_family,
    patches       = PARAMS$n_patches,
    patch_size    = PARAMS$patch_size_ha,
    env_autoscale = PARAMS$env_autoscale,
    device        = PARAMS$device,
    folder        = NULL,        # no checkpoint files ...
    checkpoints   = 10L,         # ... but a parameter snapshot every 10 epochs in model$param_history
    plot_progress = FALSE)
fit_minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

## ---- 4. save: the parameter list first, then everything else ---------------
parameters <- lapply(model$parameters, function(p) p$detach()$cpu()$clone())
torch_save(parameters, file.path(FIT_DIR, "parameters.pt"))
torch_save(model, file.path(FIT_DIR, "model.pt"))
saveRDS(model$param_history, file.path(FIT_DIR, "param_history.rds"))
loss_hist <- as.data.table(do.call(rbind, model$history))
setnames(loss_hist, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(loss_hist))])
loss_hist[, epoch := .I]
fwrite(loss_hist, file.path(FIT_DIR, "loss_history.csv"))
fwrite(loss_hist, file.path(CAL_DIR, "loss_history.csv"))
fwrite(as.data.table(model$env_scaling), file.path(CAL_DIR, "env_scaling.csv"))
## both files are release assets: the name-matched parameter list (what the
## package rebuilds from) and the whole model object (torch_load() under the
## same FINN version restores it directly)
sha <- digest::digest(file = file.path(FIT_DIR, "parameters.pt"), algo = "sha256")
writeLines(sha, file.path(CAL_DIR, "parameters_sha256.txt"))
writeLines(digest::digest(file = file.path(FIT_DIR, "model.pt"), algo = "sha256"), file.path(CAL_DIR, "model_sha256.txt"))

## ---- 5. the package must rebuild exactly this model from the list ----------
## R/definition.R is the user-facing copy of section 2; if the two drift apart
## the released parameters would silently drive a different model.
source(file.path("R", "definition.R"))
rebuilt <- build_model(variant = VARIANT, parameters = torch_load(file.path(FIT_DIR, "parameters.pt")),
                       env_scaling = as.data.frame(model$env_scaling))
model$eval(); rebuilt$eval()
chk_env <- env[siteID %in% head(unique(siteID), 20)]
chk_init <- makeInitCohorts(init[siteID %in% unique(chk_env$siteID)], Nspecies = n_species)
FINN.seed(1); s1 <- simulateForest(model,   env = chk_env, init_cohort = chk_init, patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")
FINN.seed(1); s2 <- simulateForest(rebuilt, env = chk_env, init_cohort = chk_init, patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")
max_diff <- max(abs(s1$long$site$value - s2$long$site$value), na.rm = TRUE)
if (!is.finite(max_diff) || max_diff > 1e-6) stop("build_model() does not reproduce the fitted model (max |diff| = ", max_diff, ")")

## ---- 6. fit_summary.md -----------------------------------------------------
K <- exp(as.numeric(as_array(model$reg_logK)))
writeLines(c(
  sprintf("# calibration/%s: fit summary", VARIANT), "",
  sprintf("Fitted %s by `calibration/%s/02_fit.R` in %.0f min (%.1f s per epoch, %d thread) on %s (FINN %s, torch %s, R %s).",
          format(Sys.time(), "%Y-%m-%d %H:%M"), VARIANT, fit_minutes, 60 * fit_minutes / PARAMS$epochs, PARAMS$threads, PARAMS$device,
          as.character(packageVersion("FINN")), as.character(packageVersion("torch")), paste(R.version$major, R.version$minor, sep = ".")), "",
  sprintf("Data: `%s/%s` (%d sites, %d species incl. other, 2 census steps of 10 years). Holdout for validation: `%s`.",
          PARAMS$prep_version, PARAMS$train_set, uniqueN(obs$siteID), n_species, PARAMS$holdout_set), "",
  sprintf("Processes: FINN::growth, FINN::mortality, FINN::regeneration_saturation (K init %g), FINN::competition; climate `%s` linear per species (hidden = NULL); reg_floor %g; recruits_dbh %g; %d patches of %g ha.",
          PARAMS$reg_K_init, paste(deparse(PARAMS$env_formula), collapse = ""), PARAMS$reg_floor, PARAMS$recruits_dbh, PARAMS$n_patches, PARAMS$patch_size_ha), "",
  sprintf("Fit: %d epochs, batch %d, lr %g, weights \"%s\" with weight_scale \"%s\" (baseline-scaled losses), loss families %s, env_autoscale %s, seed %d.",
          PARAMS$epochs, PARAMS$batchsize, PARAMS$lr, PARAMS$weights, PARAMS$weight_scale,
          paste(sprintf("%s = %s", names(PARAMS$loss_family), PARAMS$loss_family), collapse = ", "), PARAMS$env_autoscale, PARAMS$seed), "",
  sprintf("Final losses (baseline-scaled, 1 = no better than the mean): %s.",
          paste(sprintf("%s %.3f", names(loss_hist)[1:6], as.numeric(loss_hist[.N, 1:6])), collapse = ", ")), "",
  sprintf("Fitted recruitment caps K (stems per ha per decade): min %.1f, median %.1f, max %.1f.", min(K), median(K), max(K)), "",
  sprintf("Environment scaling stored with the model: %s.",
          paste(sprintf("%s (center %.3f, scale %.3f)", model$env_scaling$variable, model$env_scaling$center, model$env_scaling$scale), collapse = "; ")), "",
  sprintf("parameters.pt: %d parameter tensors (%s); sha256 %s. build_model() reproduces the fitted simulation (max |diff| %.2g).",
          length(parameters), paste(names(parameters), collapse = ", "), sha, max_diff), "",
  sprintf("Validation report: `Rscript calibration/%s/03_validation.R`.", VARIANT)
), file.path(CAL_DIR, "fit_summary.md"))
message(sprintf("done in %.0f min; final losses: %s", fit_minutes,
                paste(sprintf("%s=%.3f", names(loss_hist)[1:6], as.numeric(loss_hist[.N, 1:6])), collapse = " ")))
