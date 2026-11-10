#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/fit_oregon.R   (variant "oregon", v0.2 production fit)
##
## WHAT   The v0.2 production fit of the Oregon model: the v0.1 structure
##        (FINN::growth, FINN::mortality, FINN::competition,
##        FINN::regeneration_saturation, linear climate response on temp and
##        prec, 10-year step) with the settings the v0.2 audit and pilots
##        motivated:
##        1. DATA v2 (oregon-10yr-v2): the lumped `other` class dropped
##           (audit G: 1 % of the observed recruits, 27 % of the predicted ones,
##           K 1.3e5).
##        2. RECRUITMENT CAP fixed per species at 2 x the 95th percentile of
##           the observed recruits per ha per decade (floor 20), not trainable
##           (audit F2: K never bound; FINN cannot bound a custom parameter).
##        3. DISPERSION theta frozen at the v0.1 values clipped to [0.05, 2]
##           (brief trap 7, audit F3).
##        4. STAGED FIT (brief trap 9, audit E4/E5): stage A demography (states,
##           growth, mortality) with the recruitment process frozen; stage B
##           recruitment alone, everything else frozen; cosine lr in each stage
##           (audit F1: the v0.1 losses jumped in the last fifth).
##        SETTINGS TO CONFIRM FROM THE PILOTS (calibration/v02/PILOTS_v02.md):
##        if O2 (joint, capped, frozen theta) beats O3 (staged) on the holdout
##        change metrics and the bare-ground run, PARAMS$staged is set FALSE.
## DATA   <DATA_DIR>/prepared/oregon-10yr-v2/train_blocked (fit), holdout_blocked (validation only)
## RUN    Rscript calibration/v02/fit_oregon.R        (repo root; one thread; ~2.2 s/epoch on the
##        laptop -> 400 + 200 epochs ~ 25 min). NOT started by the agent (v0.2 brief).
## OUTPUT <DATA_DIR>/fits/oregon-v02/{parameters.pt, model.pt, param_history_A.rds, param_history_B.rds, loss_history.csv}
##        calibration/v02/oregon/{fit_summary.md, loss_history.csv, env_scaling.csv, parameters_sha256.txt,
##        model_sha256.txt, reg_cap.csv, recruitment_parameters.csv}
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
source(file.path("R", "definition.R"))
torch_set_num_threads(1L)
## torch >= 0.14 has the faster ignite Adam (FINN's default); the cluster's torch 0.13 falls back to optim_adam
OPTIMIZER <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
DATA_DIR <- normalizePath(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))[
  dir.exists(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia")))][1])
VARIANT <- "oregon"
FIT_DIR <- file.path(DATA_DIR, "fits", "oregon-v02"); dir.create(FIT_DIR, recursive = TRUE, showWarnings = FALSE)
CAL_DIR <- file.path("calibration", "v02", VARIANT); dir.create(CAL_DIR, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  prep_version  = "oregon-10yr-v2", train_set = "train_blocked", holdout_set = "holdout_blocked",
  patch_size_ha = 0.01681, n_patches = 4L, recruits_dbh = 12.9,
  env_formula   = ~ temp + prec, hidden = NULL,
  cap_quantile  = 0.95, cap_factor = 2, cap_floor = 20,    # K per ha per decade, fixed (audit F2)
  reg_floor     = 0.01,
  theta_source  = "oregon", theta_bounds = c(0.05, 2),      # frozen dispersion from the v0.1 fit, species 1..27 (audit F3)
  staged        = TRUE,                                     # stage A demography, stage B recruitment (audit E4/E5; pilot O3)
  epochs_A      = 400L, epochs_B = 200L,
  batchsize     = 50L, lr = 0.01, lr_scheduler = "cosine",  # cosine anneal within each stage (audit F1)
  weights_A     = c(1, 1, 1, 1, 1, 0), weights_B = c(0, 0, 0, 0, 0, 1), weight_scale = "auto",
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  seed          = 42L
)

## ---- 1. data ---------------------------------------------------------------
train_dir <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs <- fread(file.path(train_dir, "obs_dt.csv")); env <- fread(file.path(train_dir, "env_dt.csv"))
init <- fread(file.path(train_dir, "init_cohort_dt.csv"))
species <- fread(file.path(train_dir, "species_dt.csv"))[order(species)]; n_species <- nrow(species)
init_cohort <- makeInitCohorts(init, Nspecies = n_species)
cap <- obs[, .(cap = pmax(PARAMS$cap_floor, PARAMS$cap_factor * quantile(reg, PARAMS$cap_quantile, na.rm = TRUE))), by = species]
cap <- merge(species[, .(species, species_name)], cap, by = "species", all.x = TRUE)[is.na(cap), cap := PARAMS$cap_floor][order(species)]
fwrite(cap, file.path(CAL_DIR, "reg_cap.csv"))
m01 <- build_model(PARAMS$theta_source, torch_load(file.path(DATA_DIR, "fits", PARAMS$theta_source, "parameters.pt")),
                   env_scaling = as.data.frame(fread(file.path("calibration", PARAMS$theta_source, "env_scaling.csv"))))
theta0 <- pmin(pmax(as.numeric(as_array(m01$par_theta_recruits))[seq_len(n_species)], PARAMS$theta_bounds[1]), PARAMS$theta_bounds[2])
message(sprintf("train: %d sites, %d species; K cap median %.0f (max %.0f) per ha per decade; theta median %.3f",
                uniqueN(obs$siteID), n_species, median(cap$cap), max(cap$cap), median(theta0)))

## ---- 2. the model, written out ---------------------------------------------
FINN.seed(PARAMS$seed)
model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(PARAMS$env_formula, FINN::growth, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(PARAMS$env_formula, FINN::mortality, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(PARAMS$env_formula, FINN::regeneration_saturation, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       dispersion_parameter = theta0,
                                       custom_parameters = list(reg_logK = list(init = log(cap$cap), trainable = FALSE))))
model$par_theta_recruits_raw$requires_grad_(FALSE)
dummy <- data.table(siteID = 1L, year = 1L, temp = mean(env$temp), prec = mean(env$prec))
invisible(simulateForest(model, env = dummy, init_cohort = NULL, patches = 1L, patch_size = PARAMS$patch_size_ha, device = "cpu"))
model$env_scaling <- NULL
set_trainable <- function(model, pattern, on) for (nm in names(model$parameters)) if (grepl(pattern, nm)) model$parameters[[nm]]$requires_grad_(on)
fit_stage <- function(epochs, weights) fit(model, env = env, data = obs, init_cohort = init_cohort,
    epochs = epochs, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler, lr_scheduler_params = list(T_max = epochs),
    weights = weights, weight_scale = PARAMS$weight_scale, loss_family = PARAMS$loss_family,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, env_autoscale = TRUE, device = "cpu",
    optimizer = OPTIMIZER, folder = NULL, checkpoints = 10L, plot_progress = FALSE)
loss_tab <- function(stage) { lh <- as.data.table(do.call(rbind, model$history)); setnames(lh, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(lh))]); lh[, `:=`(epoch = .I, stage = stage)][] }

## ---- 3. fit ----------------------------------------------------------------
t0 <- Sys.time()
if (PARAMS$staged) {
  set_trainable(model, "^nn_regeneration|^par_regeneration|^par_theta", FALSE)
  fit_stage(PARAMS$epochs_A, PARAMS$weights_A)
  lhA <- loss_tab("A"); saveRDS(model$param_history, file.path(FIT_DIR, "param_history_A.rds"))
  torch_save(model, file.path(FIT_DIR, "model_stage_A.pt"))
  set_trainable(model, "^nn_regeneration|^par_regeneration", TRUE)                 # theta stays frozen; K is a fixed buffer
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
parameters$reg_logK <- model$reg_logK$detach()$cpu()$clone()      # the fixed cap is a buffer, not a parameter; the list carries it
torch_save(parameters, file.path(FIT_DIR, "parameters.pt"))
torch_save(model, file.path(FIT_DIR, "model.pt"))
fwrite(loss_hist, file.path(FIT_DIR, "loss_history.csv")); fwrite(loss_hist, file.path(CAL_DIR, "loss_history.csv"))
fwrite(as.data.table(model$env_scaling), file.path(CAL_DIR, "env_scaling.csv"))
writeLines(digest::digest(file = file.path(FIT_DIR, "parameters.pt"), algo = "sha256"), file.path(CAL_DIR, "parameters_sha256.txt"))
writeLines(digest::digest(file = file.path(FIT_DIR, "model.pt"), algo = "sha256"), file.path(CAL_DIR, "model_sha256.txt"))
reg_par <- merge(species, data.table(species = seq_len(n_species), K = exp(as.numeric(as_array(model$reg_logK))),
                                     theta = as.numeric(as_array(model$par_theta_recruits))), by = "species")
fwrite(reg_par, file.path(CAL_DIR, "recruitment_parameters.csv"))
rebuilt <- build_model("oregon_v02", torch_load(file.path(FIT_DIR, "parameters.pt")), env_scaling = as.data.frame(model$env_scaling))
model$eval(); rebuilt$eval()
chk_env <- env[siteID %in% head(unique(siteID), 20)]; chk_init <- makeInitCohorts(init[siteID %in% unique(chk_env$siteID)], Nspecies = n_species)
FINN.seed(1); s1 <- simulateForest(model,   env = chk_env, init_cohort = chk_init, patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")
FINN.seed(1); s2 <- simulateForest(rebuilt, env = chk_env, init_cohort = chk_init, patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")
max_diff <- max(abs(s1$long$site$value - s2$long$site$value), na.rm = TRUE)
if (!is.finite(max_diff) || max_diff > 1e-6) stop("build_model('oregon_v02') does not reproduce the fitted model (max |diff| = ", max_diff, ")")
writeLines(c(
  "# calibration/v02/oregon: fit summary", "",
  sprintf("Fitted %s by `calibration/v02/fit_oregon.R` in %.0f min (%.1f s per epoch, 1 thread) on cpu (FINN %s, torch %s, R %s).",
          format(Sys.time(), "%Y-%m-%d %H:%M"), fit_minutes, 60 * fit_minutes / (PARAMS$epochs_A + PARAMS$epochs_B),
          as.character(packageVersion("FINN")), as.character(packageVersion("torch")), paste(R.version$major, R.version$minor, sep = ".")), "",
  sprintf("Data: `%s/%s` (%d sites, %d species, no `other`). Holdout: `%s`.", PARAMS$prep_version, PARAMS$train_set, uniqueN(obs$siteID), n_species, PARAMS$holdout_set), "",
  sprintf("Processes as v0.1 (regeneration_saturation, linear climate `%s`); K FIXED per species at %g x q%.0f of the observed recruits per ha per decade (floor %g; median %.0f, max %.0f); theta FROZEN at the v0.1 values clipped to [%g, %g] (median %.3f); reg_floor %g; recruits_dbh %g; %d patches of %g ha.",
          paste(deparse(PARAMS$env_formula), collapse = ""), PARAMS$cap_factor, 100 * PARAMS$cap_quantile, PARAMS$cap_floor, median(cap$cap), max(cap$cap),
          PARAMS$theta_bounds[1], PARAMS$theta_bounds[2], median(theta0), PARAMS$reg_floor, PARAMS$recruits_dbh, PARAMS$n_patches, PARAMS$patch_size_ha), "",
  if (PARAMS$staged) sprintf("Fit: STAGED. Stage A %d epochs, weights (%s) with the recruitment process frozen; stage B %d epochs, weights (%s), recruitment alone. Batch %d, lr %g with cosine decay in each stage, weight_scale \"%s\" (the numeric weights act on the baseline-scaled losses), env_autoscale TRUE, seed %d.",
                             PARAMS$epochs_A, paste(PARAMS$weights_A, collapse = ","), PARAMS$epochs_B, paste(PARAMS$weights_B, collapse = ","), PARAMS$batchsize, PARAMS$lr, PARAMS$weight_scale, PARAMS$seed)
  else sprintf("Fit: JOINT, %d epochs, weights auto, batch %d, lr %g cosine, seed %d.", PARAMS$epochs_A + PARAMS$epochs_B, PARAMS$batchsize, PARAMS$lr, PARAMS$seed), "",
  sprintf("Final losses (baseline-scaled): %s.", paste(sprintf("%s %.3f", names(loss_hist)[1:6], as.numeric(loss_hist[.N, 1:6])), collapse = ", ")), "",
  sprintf("Environment scaling: %s.", paste(sprintf("%s (center %.3f, scale %.3f)", model$env_scaling$variable, model$env_scaling$center, model$env_scaling$scale), collapse = "; ")), "",
  sprintf("parameters.pt: %d tensors (incl. the fixed cap reg_logK), sha256 %s; model.pt sha256 %s. build_model('oregon_v02') reproduces the fitted simulation (max |diff| %.2g).",
          length(parameters), readLines(file.path(CAL_DIR, "parameters_sha256.txt")), readLines(file.path(CAL_DIR, "model_sha256.txt")), max_diff), "",
  "Validation report: `Rscript calibration/v02/validation_oregon.R`."
), file.path(CAL_DIR, "fit_summary.md"))
message(sprintf("done in %.0f min", fit_minutes))
