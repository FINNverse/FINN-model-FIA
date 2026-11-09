#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/conus/02_fit.R   (variant "conus")
##
## WHAT   Fits the conterminous-US FIA model: FINN with mechanistic growth,
##        mortality and competition, recruitment limited by the density of
##        conspecific adults (FINN::regeneration_adult, Beverton-Holt capped),
##        and one small climate network per process (hidden 20-20, SELU,
##        dropout 0.1) on yearly TerraClimate temperature, precipitation and
##        the water-deficit anomaly of the year. Annual time step; each plot
##        keeps its real remeasurement intervals (period_length). Growth and
##        mortality are scored PER TREE (tree_data + loss_aggregation = "tree"),
##        basal area, stem numbers, mean dbh and recruitment as site means.
## WHY    This is the configuration FINN-fia validated as model B
##        (fits/2026-09-10-xy-20sp-adultreg-treeloss, 2026-09): per-tree losses
##        remove the compositional confound of stand-mean rates (the synthetic
##        recovery test halved the growth and mortality parameter error), the
##        adult-density term gives recruitment a seed source so the climate
##        network no longer has to explain from climate alone where a species
##        already grows. Differences from FINN-fia, for the tool: (1) plot
##        coordinates x, y are NOT network inputs (FINN-fia found the holdout
##        skill with x, y at the training mean equal to that with real
##        coordinates, so they add nothing a user could transfer); (2) the
##        environment is passed in natural units and standardised by FINN
##        itself, so users supply natural units; (3) FINN 0.2.0 main
##        (0848b60) instead of the loss-aggregation branch. FINN-fia's later
##        multiple-shooting variants (D-F) need per-site anchors that
##        `fit(reanchor = , anchor_cohorts = )` in 0848b60 keys by global
##        timestep, which irregular census years cannot satisfy; they are not
##        reproducible under the release FINN and are not used.
## DATA   <DATA_DIR>/prepared/conus-annual-v1/train    (training: 1500 plots, 20 species)
##        <DATA_DIR>/prepared/conus-annual-v1/holdout  (validation only, 03_validation.R)
## RUN    Rscript calibration/conus/02_fit.R        (from the repository root)
##        FINN-fia measured 19.7 s per epoch on 8 laptop cores for this
##        configuration (400 epochs = 131 min); single-threaded expect several
##        hours. The cluster wrapper calibration/hpc/run_fit.sbatch passes only this path.
## OUTPUT <DATA_DIR>/fits/conus/parameters.pt   named list of parameter tensors (the release asset)
##        <DATA_DIR>/fits/conus/model.pt        the whole fitted object (validation convenience, not distributed)
##        <DATA_DIR>/fits/conus/param_history.rds, loss_history.csv
##        calibration/conus/fit_summary.md, loss_history.csv, env_scaling.csv, parameters_sha256.txt,
##        adult_parameters.csv
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
VARIANT <- "conus"
FIT_DIR <- file.path(DATA_DIR, "fits", VARIANT); dir.create(FIT_DIR, recursive = TRUE, showWarnings = FALSE)
CAL_DIR <- file.path("calibration", VARIANT)

PARAMS <- list(
  prep_version  = "conus-annual-v1",
  train_set     = "train",
  holdout_set   = "holdout",               # not read here; recorded for the validation script
  patch_size_ha = 0.01681,                 # one FIA subplot, 24 ft radius
  n_patches     = 4L,                      # the four subplots of a plot
  recruits_dbh  = 12.9,                    # a recruit enters just above the 12.7 cm tally threshold
  env_formula   = ~ tavg_c + prec_mm + cwd_z,   # yearly TerraClimate: temperature (C), precipitation (mm), deficit z-anomaly
  hidden        = c(20L, 20L),             # climate network per process, as in FINN-fia model B
  dropout       = 0.1,
  reg_K_init    = 50,                      # Beverton-Holt cap init, stems/ha/yr (FINN-fia)
  reg_floor     = 1e-3,                    # recruitment without a conspecific adult (immigration), stems/ha/yr
  ## adult threshold init: 20 cm lies between the lower quartile and the median
  ## initial dbh of most species, so x receives gradient from the start
  adult_x_init  = 20,
  ## half-saturation init: one adult on the plot is 1 / (4 x 0.01681 ha) = 14.9
  ## stems/ha, so recruitment is at half its maximum with a single seed tree
  adult_A_init  = 15,
  epochs        = 400L,                    # FINN-fia model B; its continuation to 2000 epochs did not improve the holdout
  batchsize     = 100L,                    # 15 batches per epoch
  lr            = 0.01,
  ## every term is divided by its intercept-only baseline before the terms are
  ## summed, so the six responses start at equal weight (FINN >= 0.2.0 default);
  ## a numeric vector would multiply these baseline-scaled terms
  weights       = "auto",
  weight_scale  = "auto",
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse",
                    mortality = "binomial", regeneration = "nbinom"),
  loss_aggregation = list(growth = "tree", mortality = "tree"),   # scored on the observed trees
  env_autoscale = TRUE,                    # natural units in, z-scaling learned and stored with the model
  device        = "cpu",
  threads       = 1L,                      # single-threaded: the fits share a node
  seed          = 1L
)
torch_set_num_threads(PARAMS$threads)

## ---- 1. data ---------------------------------------------------------------
train_dir <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs     <- fread(file.path(train_dir, "obs_dt.csv"))          # site x census x species; period_length sets the annual windows
env     <- fread(file.path(train_dir, "env_dt.csv"))          # site x year, natural units
init    <- fread(file.path(train_dir, "init_cohort_dt.csv"))  # trees at the first inventory
trees   <- fread(file.path(train_dir, "tree_obs_dt.csv"))     # one row per tree per interval
species <- fread(file.path(train_dir, "species_dt.csv"))[order(species)]
n_species <- nrow(species)
init_cohort <- makeInitCohorts(init, Nspecies = n_species)
message(sprintf("train: %d sites, %d species, census years %d..%d, %d observation rows; %d tree intervals (%d with growth, %d deaths)",
                uniqueN(obs$siteID), n_species, min(obs$year), max(obs$year), nrow(obs),
                nrow(trees), sum(is.finite(trees$growth)), sum(trees$died)))

## ---- 2. the model, written out ---------------------------------------------
FINN.seed(PARAMS$seed)
model <- finn(
  N_species    = n_species,
  recruits_dbh = PARAMS$recruits_dbh,
  reg_floor    = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(PARAMS$env_formula, FINN::growth,
                                       hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(PARAMS$env_formula, FINN::mortality,
                                       hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(PARAMS$env_formula, FINN::regeneration_adult,
                                       hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       custom_parameters = list(
                                         reg_logK       = rep(log(PARAMS$reg_K_init),   n_species),
                                         reg_adult_logx = rep(log(PARAMS$adult_x_init), n_species),
                                         reg_adult_logA = rep(log(PARAMS$adult_A_init), n_species)))
)

## ---- 3. fit ----------------------------------------------------------------
t0 <- Sys.time()
fit(model,
    env              = env,
    data             = obs,
    tree_data        = trees,
    init_cohort      = init_cohort,
    epochs           = PARAMS$epochs,
    batchsize        = PARAMS$batchsize,
    lr               = PARAMS$lr,
    weights          = PARAMS$weights,
    weight_scale     = PARAMS$weight_scale,
    loss_family      = PARAMS$loss_family,
    loss_aggregation = PARAMS$loss_aggregation,
    patches          = PARAMS$n_patches,
    patch_size       = PARAMS$patch_size_ha,
    env_autoscale    = PARAMS$env_autoscale,
    device           = PARAMS$device,
    folder           = NULL,     # no checkpoint files ...
    checkpoints      = 10L,      # ... but a parameter snapshot every 10 epochs in model$param_history
    plot_progress    = FALSE)
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
sha <- digest::digest(file = file.path(FIT_DIR, "parameters.pt"), algo = "sha256")
writeLines(sha, file.path(CAL_DIR, "parameters_sha256.txt"))
adult_par <- data.table(species = seq_len(n_species),
                        x_cm   = exp(as.numeric(as_array(model$reg_adult_logx))),
                        A_half = exp(as.numeric(as_array(model$reg_adult_logA))),
                        K      = exp(as.numeric(as_array(model$reg_logK))))
adult_par <- merge(species[, .(species, species_name)], adult_par, by = "species")
fwrite(adult_par, file.path(CAL_DIR, "adult_parameters.csv"))

## ---- 5. the package must rebuild exactly this model from the list ----------
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
writeLines(c(
  sprintf("# calibration/%s: fit summary", VARIANT), "",
  sprintf("Fitted %s by `calibration/%s/02_fit.R` in %.0f min (%.1f s per epoch, %d thread) on %s (FINN %s, torch %s, R %s).",
          format(Sys.time(), "%Y-%m-%d %H:%M"), VARIANT, fit_minutes, 60 * fit_minutes / PARAMS$epochs, PARAMS$threads, PARAMS$device,
          as.character(packageVersion("FINN")), as.character(packageVersion("torch")), paste(R.version$major, R.version$minor, sep = ".")), "",
  sprintf("Data: `%s/%s` (%d sites, %d species, annual step, per-site intervals, %d tree intervals). Holdout for validation: `%s`.",
          PARAMS$prep_version, PARAMS$train_set, uniqueN(obs$siteID), n_species, nrow(trees), PARAMS$holdout_set), "",
  sprintf("Processes: FINN::growth, FINN::mortality, FINN::regeneration_adult (inits K %g stems/ha/yr, adult threshold %g cm, half-saturation %g adults/ha), FINN::competition; climate network hidden = c(%s), dropout %g, on `%s`; reg_floor %g; recruits_dbh %g; %d patches of %g ha.",
          PARAMS$reg_K_init, PARAMS$adult_x_init, PARAMS$adult_A_init, paste(PARAMS$hidden, collapse = ", "), PARAMS$dropout,
          paste(deparse(PARAMS$env_formula), collapse = ""), PARAMS$reg_floor, PARAMS$recruits_dbh, PARAMS$n_patches, PARAMS$patch_size_ha), "",
  sprintf("Fit: %d epochs, batch %d, lr %g, weights \"%s\" with weight_scale \"%s\" (baseline-scaled losses), loss families %s, growth and mortality aggregated per tree, env_autoscale %s, seed %d.",
          PARAMS$epochs, PARAMS$batchsize, PARAMS$lr, PARAMS$weights, PARAMS$weight_scale,
          paste(sprintf("%s = %s", names(PARAMS$loss_family), PARAMS$loss_family), collapse = ", "), PARAMS$env_autoscale, PARAMS$seed), "",
  sprintf("Final losses (baseline-scaled, 1 = no better than the mean): %s.",
          paste(sprintf("%s %.3f", names(loss_hist)[1:6], as.numeric(loss_hist[.N, 1:6])), collapse = ", ")), "",
  sprintf("Fitted adult threshold x: min %.1f, median %.1f, max %.1f cm. Half-saturation: min %.1f, median %.1f, max %.1f adults/ha. Caps K: median %.1f stems/ha/yr (adult_parameters.csv).",
          min(adult_par$x_cm), median(adult_par$x_cm), max(adult_par$x_cm),
          min(adult_par$A_half), median(adult_par$A_half), max(adult_par$A_half), median(adult_par$K)), "",
  sprintf("Environment scaling stored with the model: %s.",
          paste(sprintf("%s (center %.3f, scale %.3f)", model$env_scaling$variable, model$env_scaling$center, model$env_scaling$scale), collapse = "; ")), "",
  sprintf("parameters.pt: %d parameter tensors (%s); sha256 %s. build_model() reproduces the fitted simulation (max |diff| %.2g).",
          length(parameters), paste(names(parameters), collapse = ", "), sha, max_diff), "",
  sprintf("Validation report: `Rscript calibration/%s/03_validation.R`.", VARIANT)
), file.path(CAL_DIR, "fit_summary.md"))
message(sprintf("done in %.0f min; final losses: %s", fit_minutes,
                paste(sprintf("%s=%.3f", names(loss_hist)[1:6], as.numeric(loss_hist[.N, 1:6])), collapse = " ")))
