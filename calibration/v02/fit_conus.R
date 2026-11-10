#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/fit_conus.R   (variant "conus", v0.2 production fit)
##
## WHAT   The v0.2 production fit of the conterminous-US model: the v0.1 model
##        structure (FINN::growth, FINN::mortality, FINN::competition,
##        FINN::regeneration_adult, one 20-20 climate network per process on
##        yearly tavg_c, prec_mm, cwd_z, per-tree growth and mortality losses)
##        with the settings the v0.2 audit and pilots motivated:
##        1. DATA v2 (conus-annual-v2): the complete ingrowth target (audit B).
##        2. RECRUITMENT CAP fixed per species at 2 x the 95th percentile of
##           the observed annual recruits per ha (floor 2), not trainable
##           (audit F2; FINN cannot bound a custom parameter).
##        3. DISPERSION theta per species frozen at the v0.1 values clipped to
##           [0.05, 2] (brief trap 7, audit F3).
##        4. STAGED FIT (brief trap 9, audit E1/E3): stage A demography (per-tree
##           growth and mortality, states) with the recruitment process frozen;
##           stage B recruitment alone, everything else frozen. Cosine lr in
##           each stage (audit F1).
##        SETTINGS TO CONFIRM FROM THE PILOTS (calibration/v02/PILOTS_v02.md):
##        the stage split; if C1 (joint, capped) beats C2 (staged) on the
##        holdout change metrics and AUC, PARAMS$staged is set FALSE and the
##        fit is one joint run of PARAMS$epochs_A + PARAMS$epochs_B epochs.
## DATA   <DATA_DIR>/prepared/conus-annual-v2/train (fit), holdout (validation only)
## RUN    Rscript calibration/v02/fit_conus.R        (repo root; one thread; ~35 s/epoch on
##        the laptop -> 400 + 150 epochs ~ 5.5 h). NOT started by the agent (v0.2 brief).
## OUTPUT <DATA_DIR>/fits/conus-v02/{parameters.pt, model.pt, param_history_A.rds, param_history_B.rds, loss_history.csv}
##        calibration/v02/conus/{fit_summary.md, loss_history.csv, env_scaling.csv, parameters_sha256.txt,
##        model_sha256.txt, adult_parameters.csv, reg_cap.csv}
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
source(file.path("R", "definition.R"))
torch_set_num_threads(1L)
## torch >= 0.14 has the faster ignite Adam (FINN's default); the cluster's torch 0.13 falls back to optim_adam
OPTIMIZER <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
DATA_DIR <- normalizePath(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))[
  dir.exists(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia")))][1])
VARIANT <- "conus"
FIT_DIR <- file.path(DATA_DIR, "fits", "conus-v02"); dir.create(FIT_DIR, recursive = TRUE, showWarnings = FALSE)
CAL_DIR <- file.path("calibration", "v02", VARIANT); dir.create(CAL_DIR, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  prep_version  = "conus-annual-v2", train_set = "train", holdout_set = "holdout",
  patch_size_ha = 0.01681, n_patches = 4L, recruits_dbh = 12.9,
  env_formula   = ~ tavg_c + prec_mm + cwd_z, hidden = c(20L, 20L), dropout = 0.1,
  cap_quantile  = 0.95, cap_factor = 2, cap_floor = 2,     # K per ha per year, fixed (audit F2)
  reg_floor     = 1e-3, adult_x_init = 20, adult_A_init = 15,
  theta_source  = "conus", theta_bounds = c(0.05, 2),       # frozen dispersion from the v0.1 fit (audit F3)
  staged        = TRUE,                                     # stage A demography, stage B recruitment (audit E1/E3; pilot C2)
  epochs_A      = 400L, epochs_B = 150L,                    # v0.1 length for the demography, a recruitment-only tail
  batchsize     = 100L, lr = 0.01, lr_scheduler = "cosine", # cosine anneal within each stage (audit F1)
  weights_A     = c(1, 1, 1, 1, 1, 0), weights_B = c(0, 0, 0, 0, 0, 1), weight_scale = "auto",
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  loss_aggregation = list(growth = "tree", mortality = "tree"),
  seed          = 1L
)

## ---- 1. data ---------------------------------------------------------------
train_dir <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs <- fread(file.path(train_dir, "obs_dt.csv")); env <- fread(file.path(train_dir, "env_dt.csv"))
init <- fread(file.path(train_dir, "init_cohort_dt.csv")); trees <- fread(file.path(train_dir, "tree_obs_dt.csv"))
species <- fread(file.path(train_dir, "species_dt.csv"))[order(species)]; n_species <- nrow(species)
init_cohort <- makeInitCohorts(init, Nspecies = n_species)
## the fixed cap: 2 x q95 of the observed annual recruits per ha of each species (floor 2)
obs_annual <- copy(obs)[, reg_yr := reg / period_length]
cap <- obs_annual[, .(cap = pmax(PARAMS$cap_floor, PARAMS$cap_factor * quantile(reg_yr, PARAMS$cap_quantile, na.rm = TRUE))), by = species]
cap <- merge(species[, .(species, species_name)], cap, by = "species", all.x = TRUE)[is.na(cap), cap := PARAMS$cap_floor][order(species)]
fwrite(cap, file.path(CAL_DIR, "reg_cap.csv"))
## the frozen dispersion: the v0.1 production fit's theta per species, clipped
m01 <- build_model(PARAMS$theta_source, torch_load(file.path(DATA_DIR, "fits", PARAMS$theta_source, "parameters.pt")),
                   env_scaling = as.data.frame(fread(file.path("calibration", PARAMS$theta_source, "env_scaling.csv"))))
theta0 <- pmin(pmax(as.numeric(as_array(m01$par_theta_recruits)), PARAMS$theta_bounds[1]), PARAMS$theta_bounds[2])
message(sprintf("train: %d sites, %d species, %d tree intervals; K cap median %.1f (max %.1f) per ha per year; theta median %.3f",
                uniqueN(obs$siteID), n_species, nrow(trees), median(cap$cap), max(cap$cap), median(theta0)))

## ---- 2. the model, written out ---------------------------------------------
FINN.seed(PARAMS$seed)
model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(PARAMS$env_formula, FINN::growth, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(PARAMS$env_formula, FINN::mortality, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(PARAMS$env_formula, FINN::regeneration_adult, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       dispersion_parameter = theta0,
                                       custom_parameters = list(reg_logK = list(init = log(cap$cap), trainable = FALSE),
                                                                reg_adult_logx = rep(log(PARAMS$adult_x_init), n_species),
                                                                reg_adult_logA = rep(log(PARAMS$adult_A_init), n_species))))
model$par_theta_recruits_raw$requires_grad_(FALSE)
## the climate layers are created lazily; instantiate them so a stage can freeze them
dummy <- data.table(siteID = 1L, year = 1L, tavg_c = mean(env$tavg_c), prec_mm = mean(env$prec_mm), cwd_z = 0)
invisible(simulateForest(model, env = dummy, init_cohort = NULL, patches = 1L, patch_size = PARAMS$patch_size_ha, device = "cpu"))
model$env_scaling <- NULL
set_trainable <- function(model, pattern, on) for (nm in names(model$parameters)) if (grepl(pattern, nm)) model$parameters[[nm]]$requires_grad_(on)
fit_stage <- function(epochs, weights) fit(model, env = env, data = obs, tree_data = trees, init_cohort = init_cohort,
    epochs = epochs, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler, lr_scheduler_params = list(T_max = epochs),
    weights = weights, weight_scale = PARAMS$weight_scale, loss_family = PARAMS$loss_family, loss_aggregation = PARAMS$loss_aggregation,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, env_autoscale = TRUE, device = "cpu",
    optimizer = OPTIMIZER, folder = NULL, checkpoints = 10L, plot_progress = FALSE)
loss_tab <- function(stage) { lh <- as.data.table(do.call(rbind, model$history)); setnames(lh, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(lh))]); lh[, `:=`(epoch = .I, stage = stage)][] }

## ---- 3. fit ----------------------------------------------------------------
t0 <- Sys.time()
if (PARAMS$staged) {
  set_trainable(model, "^nn_regeneration|^par_regeneration|^par_theta|^reg_", FALSE)
  fit_stage(PARAMS$epochs_A, PARAMS$weights_A)
  lhA <- loss_tab("A"); saveRDS(model$param_history, file.path(FIT_DIR, "param_history_A.rds"))
  torch_save(model, file.path(FIT_DIR, "model_stage_A.pt"))
  set_trainable(model, "^nn_regeneration|^par_regeneration|^reg_adult", TRUE)       # theta and K stay frozen
  set_trainable(model, "^nn_growth|^nn_mortality|^par_growth|^par_mortality|^par_competition", FALSE)
  fit_stage(PARAMS$epochs_B, PARAMS$weights_B)
  lhB <- loss_tab("B"); saveRDS(model$param_history, file.path(FIT_DIR, "param_history_B.rds"))
  loss_hist <- rbind(lhA, lhB)
  set_trainable(model, ".", TRUE); model$par_theta_recruits_raw$requires_grad_(FALSE)
} else {
  fit_stage(PARAMS$epochs_A + PARAMS$epochs_B, "auto")
  loss_hist <- loss_tab("joint"); saveRDS(model$param_history, file.path(FIT_DIR, "param_history_A.rds"))
}
fit_minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

## ---- 4. save ---------------------------------------------------------------
parameters <- lapply(model$parameters, function(p) p$detach()$cpu()$clone())
torch_save(parameters, file.path(FIT_DIR, "parameters.pt"))
torch_save(model, file.path(FIT_DIR, "model.pt"))
fwrite(loss_hist, file.path(FIT_DIR, "loss_history.csv")); fwrite(loss_hist, file.path(CAL_DIR, "loss_history.csv"))
fwrite(as.data.table(model$env_scaling), file.path(CAL_DIR, "env_scaling.csv"))
writeLines(digest::digest(file = file.path(FIT_DIR, "parameters.pt"), algo = "sha256"), file.path(CAL_DIR, "parameters_sha256.txt"))
writeLines(digest::digest(file = file.path(FIT_DIR, "model.pt"), algo = "sha256"), file.path(CAL_DIR, "model_sha256.txt"))
adult_par <- merge(species, data.table(species = seq_len(n_species), x_cm = exp(as.numeric(as_array(model$reg_adult_logx))),
                                       A_half = exp(as.numeric(as_array(model$reg_adult_logA))), K = exp(as.numeric(as_array(model$reg_logK))),
                                       theta = as.numeric(as_array(model$par_theta_recruits))), by = "species")
fwrite(adult_par, file.path(CAL_DIR, "adult_parameters.csv"))
## the package must rebuild exactly this model (R/definition.R variant "conus_v02")
rebuilt <- build_model("conus_v02", torch_load(file.path(FIT_DIR, "parameters.pt")), env_scaling = as.data.frame(model$env_scaling))
model$eval(); rebuilt$eval()
chk_env <- env[siteID %in% head(unique(siteID), 20)]; chk_init <- makeInitCohorts(init[siteID %in% unique(chk_env$siteID)], Nspecies = n_species)
FINN.seed(1); s1 <- simulateForest(model,   env = chk_env, init_cohort = chk_init, patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")
FINN.seed(1); s2 <- simulateForest(rebuilt, env = chk_env, init_cohort = chk_init, patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")
max_diff <- max(abs(s1$long$site$value - s2$long$site$value), na.rm = TRUE)
if (!is.finite(max_diff) || max_diff > 1e-6) stop("build_model('conus_v02') does not reproduce the fitted model (max |diff| = ", max_diff, ")")
writeLines(c(
  "# calibration/v02/conus: fit summary", "",
  sprintf("Fitted %s by `calibration/v02/fit_conus.R` in %.0f min (%.1f s per epoch, 1 thread) on cpu (FINN %s, torch %s, R %s).",
          format(Sys.time(), "%Y-%m-%d %H:%M"), fit_minutes, 60 * fit_minutes / (PARAMS$epochs_A + PARAMS$epochs_B),
          as.character(packageVersion("FINN")), as.character(packageVersion("torch")), paste(R.version$major, R.version$minor, sep = ".")), "",
  sprintf("Data: `%s/%s` (%d sites, %d species, %d tree intervals). Holdout: `%s`.", PARAMS$prep_version, PARAMS$train_set, uniqueN(obs$siteID), n_species, nrow(trees), PARAMS$holdout_set), "",
  sprintf("Processes as v0.1 (regeneration_adult, climate networks c(%s) on `%s`); K FIXED per species at %g x q%.0f of the observed annual recruits per ha (floor %g; median %.1f, max %.1f); theta FROZEN at the v0.1 values clipped to [%g, %g] (median %.3f); reg_floor %g; recruits_dbh %g; %d patches of %g ha.",
          paste(PARAMS$hidden, collapse = ", "), paste(deparse(PARAMS$env_formula), collapse = ""), PARAMS$cap_factor, 100 * PARAMS$cap_quantile, PARAMS$cap_floor,
          median(cap$cap), max(cap$cap), PARAMS$theta_bounds[1], PARAMS$theta_bounds[2], median(theta0), PARAMS$reg_floor, PARAMS$recruits_dbh, PARAMS$n_patches, PARAMS$patch_size_ha), "",
  if (PARAMS$staged) sprintf("Fit: STAGED. Stage A %d epochs, weights (%s) with the recruitment process frozen; stage B %d epochs, weights (%s), recruitment alone. Batch %d, lr %g with cosine decay in each stage, weight_scale \"%s\" (the numeric weights act on the baseline-scaled losses), growth and mortality per tree, env_autoscale TRUE, seed %d.",
                             PARAMS$epochs_A, paste(PARAMS$weights_A, collapse = ","), PARAMS$epochs_B, paste(PARAMS$weights_B, collapse = ","), PARAMS$batchsize, PARAMS$lr, PARAMS$weight_scale, PARAMS$seed)
  else sprintf("Fit: JOINT, %d epochs, weights auto, batch %d, lr %g cosine, seed %d.", PARAMS$epochs_A + PARAMS$epochs_B, PARAMS$batchsize, PARAMS$lr, PARAMS$seed), "",
  sprintf("Final losses (baseline-scaled): %s.", paste(sprintf("%s %.3f", names(loss_hist)[1:6], as.numeric(loss_hist[.N, 1:6])), collapse = ", ")), "",
  sprintf("Adult threshold x: median %.1f cm; half-saturation: median %.1f adults/ha (adult_parameters.csv).", median(adult_par$x_cm), median(adult_par$A_half)), "",
  sprintf("Environment scaling: %s.", paste(sprintf("%s (center %.3f, scale %.3f)", model$env_scaling$variable, model$env_scaling$center, model$env_scaling$scale), collapse = "; ")), "",
  sprintf("parameters.pt: %d tensors, sha256 %s; model.pt sha256 %s. build_model('conus_v02') reproduces the fitted simulation (max |diff| %.2g).",
          length(parameters), readLines(file.path(CAL_DIR, "parameters_sha256.txt")), readLines(file.path(CAL_DIR, "model_sha256.txt")), max_diff), "",
  "Validation report: `Rscript calibration/v02/validation_conus.R`."
), file.path(CAL_DIR, "fit_summary.md"))
message(sprintf("done in %.0f min", fit_minutes))
