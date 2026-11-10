#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-zoo | FINN-model-FIA/calibration/v02/recipe/oregon/fit_stageB.R
##
## WHAT   Stage B of the zoo v0.2 recipe for the Oregon variant (one simulation
##        step = one 10-year remeasurement interval): the regeneration process
##        (FINN::regeneration_saturation: light response x climate, Beverton-Holt
##        cap K per species) is fitted with competition, growth and mortality
##        FROZEN at stage A (fit_stageA.R) and the recruitment dispersion theta
##        frozen at its maximum-likelihood value. The result is the simulation
##        model of the zoo release (model.pt + parameters.pt).
##        Copied from FINN-model-IFN-Spain/calibration/v02/recipe/fit_stageB.R
##        (cb41778); the per-model block of PARAMS and the data handling differ.
## WHY    The recipe (FINN-zoo-data/RECIPE_v02.md): recruitment is identified by
##        the recruitment observations alone; fitting it jointly with the
##        demography lets the state losses pull weakly identified parameters (BWI;
##        in the v0.1 Oregon fit the lumped class and the firs took the recruits).
##        - regeneration env response: linear in the z-scaled climate (intercept
##          warm-started at the species' observed mean recruitment per step);
##        - theta: per species, the intercept-only ML estimate on the observed
##          recruit counts, clamped to [theta_min, theta_max] and NOT trained (a
##          free theta collapsed and recruitment exploded in the BWI fits);
##        - the recruitment observation enters the negative binomial as whole
##          tallied recruits (recruit_obs_weight = the area one recruit stands for:
##          the plot, 4 subplots of 0.01681 ha);
##        - losses: ALL terms (dbh, ba, trees, growth, mortality at the site level,
##          regeneration), weights "auto" -- the stage-B loss rule of the recipe;
##        - the frozen networks run without dropout (fit() puts the model in train
##          mode, which would jitter the frozen demography);
##        - cosine lr, plain Adam, one thread, light_mask TRUE, default floor.
##        DEVIATIONS from the reference (see fit_stageA.R): the recruits are drawn
##        with FINN's "refined" sampler (4 x 0.0168-ha patches), one step is a
##        decade (K in recruits per ha per STEP), no tree table (site-level growth
##        and mortality losses). K: trainable as the recipe says (the pilots' fixed
##        data-derived cap is not used); the fit summary reports the fitted K
##        against the observed recruitment quantiles so a runaway is visible (the
##        v0.1 fit's K reached 131,000 per ha per decade).
## DATA   ../../FINN-zoo-data/fia/prepared/oregon-10yr-v2/train_blocked  and
##        <FIT_DIR>/model_stageA.pt (written by fit_stageA.R)
## NEEDS  FINN 0.2.0 (commit 0848b60), torch, MASS (theta.ml)
## RUN    Rscript calibration/v02/recipe/oregon/fit_stageB.R      (from the repo root, after fit_stageA.R)
## OUTPUT <FIT_DIR>/model.pt            the simulation model (stage-A demography + stage-B regeneration)
##        <FIT_DIR>/parameters.pt       named list: state_dict, env_scaling, species, settings
##        <FIT_DIR>/model_init_stageB.pt, loss_history_stageB.csv, regeneration_parameters.csv,
##        <FIT_DIR>/convergence_stageB.csv, fit_summary_stageB.md
##        FIT_DIR = ../../FINN-zoo-data/fia/fits/oregon-10yr-v2-recipe (outside git)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
if (dir.exists("~/R/lib-finn-0848b60")) .libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({
  library(FINN)
  library(data.table)
  library(torch)
})
stopifnot(packageVersion("FINN") == "0.2.0")

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "fia")
FIT_DIR  <- file.path(DATA_DIR, "fits", "oregon-10yr-v2-recipe")

PARAMS <- list(
  ## ---- per model (must equal fit_stageA.R where both have the entry) ----
  prep_version  = "oregon-10yr-v2",
  train_set     = "train_blocked",
  holdout_set   = "holdout_blocked",       # not read here; recorded for validation.R
  env_vars      = c("temp", "prec"),       # as stage A
  patch_size_ha = 0.01681,                 # as stage A
  n_patches     = 4L,
  recruits_dbh  = 12.9,
  step_years    = 10L,                     # one simulation step = one remeasurement interval; K, theta and the rates are per step
  recruit_obs_weight_ha = 4 * 0.01681,     # the area one observed recruit stands for: the plot (four subplots), so reg x w is the tallied count and the NB scores counts
  reg_K_init    = 50,                      # Beverton-Holt cap start, recruits/ha/step (between the 95th percentile of all cells, 15, and the 99th, 89; the recruiting cells' median is 30)
  reg_sampler   = "refined",               # as stage A: exact Poisson draw of the per-patch recruits below a mean of 2
  ## ---- the recipe ----
  reg_floor     = 1e-3,                    # FINN's additive recruitment floor inside regeneration_saturation (stems/ha/step)
  theta_min     = 0.01,                    # the frozen dispersion never below this (BWI rule: theta >= 0.01)
  theta_max     = 100,                     # a species whose counts show no overdispersion (or no recruit at all) gets a near-Poisson theta
  sen_max       = 5,                       # as stage A: U-shape cap
  growth_scale  = 0.01,                    # as stage A: growth enters mortality as growth / 0.01
  frozen_dropout = 0,                      # dropout of the frozen growth and mortality networks during this fit: off, as in simulateForest()
  light_mask    = TRUE,
  n_quantiles   = 10L,
  epochs        = 150L,                    # the regeneration block settles faster than the demography (BWI: 100)
  batchsize     = 100L,
  lr            = 0.005,
  lr_scheduler  = "cosine",
  weights       = "auto",                  # every term scaled by its intercept-only baseline
  weight_scale  = "auto",
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse",
                    mortality = "binomial", regeneration = "nbinom"),   # ALL terms (the stage-B loss rule, see the header)
  loss_aggregation = NULL,                 # no tree table in this variant: growth and mortality at the site x species level
  checkpoints   = 10L,
  frozen_tol    = 1e-5,                    # a frozen block may not move by more than this
  threads       = 1L,
  device        = "cpu",
  seed          = 1L
)
torch_set_num_threads(PARAMS$threads)
stopifnot(identical(formals(FINN::mortality)$growth_scale, PARAMS$growth_scale))
PS <- PARAMS$patch_size_ha
w_obs <- PARAMS$recruit_obs_weight_ha

## ---- 1. data -----------------------------------------------------------------
train_dir <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs   <- fread(file.path(train_dir, "obs_dt.csv"))
env   <- fread(file.path(train_dir, "env_dt.csv"))[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init  <- fread(file.path(train_dir, "init_cohort_dt.csv"))
species_dt <- fread(file.path(train_dir, "species_dt.csv"))
n_species <- nrow(species_dt)
if (!"period_length" %in% names(obs)) obs[, period_length := 1L]   # one step per observation
init_cohort <- makeInitCohorts(init, Nspecies = n_species)
env_formula <- as.formula(paste("~", paste(PARAMS$env_vars, collapse = " + ")))
stopifnot(max(abs(obs$reg * w_obs - round(obs$reg * w_obs))) < 0.01)   # whole recruits in the observation space
iv <- unique(obs[, .(siteID, y_end = year, y0 = year - period_length)])[order(siteID, y_end)]

## ---- 2. stage A, frozen ----------------------------------------------------------
stageA <- torch_load(file.path(FIT_DIR, "model_stageA.pt"))
nn_pars <- function(nn) lapply(nn$parameters, function(p) as_array(p$detach()$cpu()))
nn_maxdiff <- function(a, b) { stopifnot(length(a) == length(b)); max(mapply(function(x, y) max(abs(x - y)), a, b)) }
for (nn in list(stageA$nn_growth, stageA$nn_mortality)) for (m in nn$modules) if (inherits(m, "nn_dropout")) m$p <- PARAMS$frozen_dropout
drop_p <- unlist(lapply(list(stageA$nn_growth, stageA$nn_mortality), function(nn)
  vapply(Filter(function(m) inherits(m, "nn_dropout"), nn$modules), function(m) m$p, numeric(1))))
stopifnot(length(drop_p) > 0, all(drop_p == PARAMS$frozen_dropout))
frozen <- list(competition = as.matrix(stageA$par_competition$detach()$cpu()), growth = as.matrix(stageA$par_growth$detach()$cpu()),
               mortality = as.matrix(stageA$par_mortality$detach()$cpu()), sen_raw = as.numeric(stageA$mort_sen_raw$detach()$cpu()),
               growth_env = nn_pars(stageA$nn_growth), mortality_env = nn_pars(stageA$nn_mortality))
scaling_A <- as.data.table(stageA$env_scaling)[variable %in% PARAMS$env_vars]
## the stage-A mortality kernel (identical to fit_stageA.R)
mortality_ushape <- function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = FALSE, growth = NULL, growth_scale = 0.01) {
  if (is.null(growth)) growth <- self$g
  d   <- dbh / 100
  sen <- (self$mort_sen_max * torch::torch_sigmoid(self$mort_sen_raw))[species]
  torch::torch_sigmoid(pred + growth / growth_scale * parMort[, 3][species] + light * parMort[, 1][species] +
                       parMort[, 2][species] * d + sen * d * d)
}

## ---- 3. the recruitment dispersion: intercept-only ML per species, frozen ----------
## theta of a negative binomial fitted to the observed recruit counts (reg x w) of
## every training plot x step with the species' mean as the only predictor; a
## species whose counts are all zero or show no overdispersion gets theta_max.
theta_ml <- vapply(seq_len(n_species), function(s) {
  y <- obs[species == s, round(reg * w_obs)]
  if (sum(y) == 0) return(PARAMS$theta_max)
  th <- suppressWarnings(tryCatch(as.numeric(MASS::theta.ml(y, mu = mean(y), limit = 100)), error = function(e) NA_real_))
  if (!is.finite(th)) PARAMS$theta_max else min(max(th, PARAMS$theta_min), PARAMS$theta_max)
}, numeric(1))
## warm start of the regeneration env layer: intercept = log of the species' mean
## recruitment per ha per step (the climate columns are z-scaled, so 0 slopes = the mean)
mean_annual_reg <- obs[, .(r = sum(reg) / sum(period_length)), by = species][order(species), r]
W0 <- cbind(log(mean_annual_reg + PARAMS$reg_floor), matrix(0, n_species, length(PARAMS$env_vars)))
## the observed recruitment quantiles the fitted K is judged against (per ha per step)
reg_q <- obs[, .(q95_cells = quantile(reg / period_length, 0.95), q90_recruiting = quantile((reg / period_length)[reg > 0], 0.9, na.rm = TRUE), max = max(reg / period_length)), by = species][order(species)]

## ---- 4. the stage-B model, written out ---------------------------------------------
FINN.seed(PARAMS$seed)
model <- finn(
  N_species    = n_species,
  recruits_dbh = PARAMS$recruits_dbh,
  reg_floor    = PARAMS$reg_floor,
  recruit_obs_weight = w_obs,              # the NB scores whole recruits (predicted density x w vs observed x w)
  competition_process  = createProcess(~0, FINN::competition, initSpecies = frozen$competition,
                                       optimizeSpecies = FALSE, n_quantiles = PARAMS$n_quantiles, light_mask = PARAMS$light_mask),
  growth_process       = createProcess(env_formula, FINN::growth, initSpecies = frozen$growth,
                                       NN = stageA$nn_growth, optimizeSpecies = FALSE, optimizeEnv = FALSE),
  mortality_process    = createProcess(env_formula, mortality_ushape, initSpecies = frozen$mortality,
                                       NN = stageA$nn_mortality, optimizeSpecies = FALSE, optimizeEnv = FALSE,
                                       custom_parameters = list(mort_sen_raw = list(init = frozen$sen_raw, trainable = FALSE),
                                                                mort_sen_max = list(init = PARAMS$sen_max, trainable = FALSE))),
  regeneration_process = createProcess(env_formula, FINN::regeneration_saturation,
                                       initEnv = list(W0), dispersion_parameter = theta_ml,
                                       custom_parameters = list(reg_logK = list(init = rep(log(PARAMS$reg_K_init), n_species), trainable = TRUE)),
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       regeneration_sampler = PARAMS$reg_sampler)
)
model$period_mode <- "per_site"
stopifnot(identical(model$regeneration_sampler, PARAMS$reg_sampler))
## the frozen networks were trained on stage A's z-scaled climate: give the model that
## scaling now (fit() recomputes it from the same training env; checked below), so the
## warm-up simulation and the start-of-stage-B recruitment check see the right inputs
model$env_scaling <- as.data.frame(scaling_A)
## FINN builds the env networks lazily; a two-site, two-step simulation builds them
## (so the passed modules become registered submodules that torch_save keeps), then
## the frozen networks and theta are taken out of the gradient explicitly (a
## simulation re-enables requires_grad on every parameter when it ends)
e_build <- env[siteID %in% sort(unique(siteID))[1:2] & year <= 2L]
e_build[, siteID := match(siteID, sort(unique(siteID)))]
invisible(simulateForest(model, env = e_build, init_cohort = NULL, patches = 1L, patch_size = PS, device = "cpu"))
for (nn in list(model$nn_growth, model$nn_mortality)) for (p in nn$parameters) invisible(p$requires_grad_(FALSE))
invisible(model$par_theta_recruits_raw$requires_grad_(FALSE))
stopifnot(nn_maxdiff(nn_pars(model$nn_growth), frozen$growth_env) < 1e-7, nn_maxdiff(nn_pars(model$nn_mortality), frozen$mortality_env) < 1e-7,
          !any(sapply(c(model$nn_growth$parameters, model$nn_mortality$parameters), function(p) p$requires_grad)),
          max(abs(as.numeric(model$par_theta_recruits$detach()) - theta_ml)) < 1e-4 * max(theta_ml),
          max(abs(as.matrix(model$par_mortality$detach()$cpu()) - frozen$mortality)) < 1e-4,
          max(abs(as_array(model$nn_regeneration$parameters[[1]]$detach()) - W0)) < 1e-6)
torch_save(model, file.path(FIT_DIR, "model_init_stageB.pt"))

## ---- 5. fit -------------------------------------------------------------------------
FINN.seed(PARAMS$seed)
t0 <- Sys.time()
fit(model,
    env           = env,
    data          = obs,
    init_cohort   = init_cohort,
    epochs        = PARAMS$epochs,
    batchsize     = PARAMS$batchsize,
    lr            = PARAMS$lr,
    lr_scheduler  = PARAMS$lr_scheduler,
    optimizer     = torch::optim_adam,
    weights       = PARAMS$weights,
    weight_scale  = PARAMS$weight_scale,
    loss_family   = PARAMS$loss_family,
    loss_aggregation = PARAMS$loss_aggregation,
    patches       = PARAMS$n_patches,
    patch_size    = PS,
    device        = PARAMS$device,
    folder        = NULL,
    checkpoints   = PARAMS$checkpoints,
    plot_progress = FALSE)
fit_minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

## save first, summarise later
model$eval()
torch_save(model, file.path(FIT_DIR, "model.pt"))
parameters <- list(
  state_dict  = model$state_dict(),                 # every parameter and buffer tensor, by name (networks included)
  env_scaling = model$env_scaling,                  # the z-scaling constants of the climate columns
  species     = species_dt,
  settings    = list(prep_version = PARAMS$prep_version, patch_size_ha = PS, n_patches = PARAMS$n_patches, step_years = PARAMS$step_years,
                     recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor, recruit_obs_weight = w_obs,
                     regeneration_sampler = PARAMS$reg_sampler,
                     env_formula = deparse(env_formula), sen_max = PARAMS$sen_max, growth_scale = PARAMS$growth_scale,
                     theta = theta_ml, epochs_stageB = PARAMS$epochs, recipe = "zoo v0.2 (stage A + B)",
                     finn_version = as.character(packageVersion("FINN")), torch_version = as.character(packageVersion("torch")))
)
torch_save(parameters, file.path(FIT_DIR, "parameters.pt"))
loss_hist <- as.data.table(do.call(rbind, model$history))
setnames(loss_hist, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(loss_hist))])
loss_hist[, epoch := .I]
fwrite(loss_hist, file.path(FIT_DIR, "loss_history_stageB.csv"))

## ---- 6. checks and what was fitted ------------------------------------------------------
## the frozen blocks must still be stage A's, the climate scaling the same
dev <- c(growth_env    = nn_maxdiff(nn_pars(model$nn_growth), frozen$growth_env),
         mortality_env = nn_maxdiff(nn_pars(model$nn_mortality), frozen$mortality_env),
         growth_par    = max(abs(as.matrix(model$par_growth$detach()$cpu()) - frozen$growth)),
         mortality_par = max(abs(as.matrix(model$par_mortality$detach()$cpu()) - frozen$mortality)),
         competition   = max(abs(as.matrix(model$par_competition$detach()$cpu()) - frozen$competition)),
         sen           = max(abs(as.numeric(model$mort_sen_raw$detach()$cpu()) - frozen$sen_raw)),
         theta_rel     = max(abs(as.numeric(model$par_theta_recruits$detach()) - theta_ml) / theta_ml))
scaling_B <- as.data.table(model$env_scaling)[variable %in% PARAMS$env_vars]
dev["env_scaling"] <- max(abs(merge(scaling_A, scaling_B, by = "variable")[, c(center.x - center.y, scale.x - scale.y)]))
message(sprintf("max |change| of the frozen blocks: %s", paste(sprintf("%s %.2g", names(dev), dev), collapse = ", ")))
stopifnot(all(dev < PARAMS$frozen_tol))
## the saved simulation model survives torch_load with its networks
chk <- torch_load(file.path(FIT_DIR, "model.pt"))
stopifnot(nn_maxdiff(nn_pars(chk$nn_growth), frozen$growth_env) < 1e-7, nn_maxdiff(nn_pars(chk$nn_mortality), frozen$mortality_env) < 1e-7,
          setequal(names(parameters$state_dict), names(chk$state_dict())), identical(chk$regeneration_sampler, PARAMS$reg_sampler))
rm(chk)
## the fitted regeneration per species
W3 <- as_array(model$nn_regeneration$parameters[[1]]$detach()$cpu())
REG <- data.table(species = seq_len(n_species),
                  mean_step_reg_obs_ha = mean_annual_reg,
                  theta = theta_ml,
                  K_init = PARAMS$reg_K_init, K_fit = exp(as.numeric(model$reg_logK$detach()$cpu())),
                  obs_q95_cells = reg_q$q95_cells, obs_q90_recruiting = reg_q$q90_recruiting, obs_max = reg_q$max,
                  light = as.numeric(as.matrix(model$par_regeneration$detach()$cpu())[, 1]),
                  `colnames<-`(as.data.table(W3), paste0("env_", c("intercept", PARAMS$env_vars))))
REG <- merge(species_dt, REG, by = "species")
fwrite(REG, file.path(FIT_DIR, "regeneration_parameters.csv"))
## predicted / observed recruitment on the training plots over fit()'s windows
## (output step k is driven by env row k; the observation at y_end is row y_end)
rec_ratio <- function(m) {
  m$eval(); FINN.seed(PARAMS$seed)
  S <- as.data.table(simulateForest(m, env = env, init_cohort = init_cohort, patches = PARAMS$n_patches, patch_size = PS, device = "cpu")$wide$site)
  for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
  W <- iv[, .(year = seq.int(y0 + 1L, y_end)), by = .(siteID, y_end)]
  sum(merge(W, S[, .(siteID, year, r_mean_ha)], by = c("siteID", "year"), allow.cartesian = TRUE)$r_mean_ha) / sum(pmax(obs$reg, 0))
}
ratio_start <- rec_ratio(torch_load(file.path(FIT_DIR, "model_init_stageB.pt")))
ratio_end   <- rec_ratio(model)
ph <- model$param_history
conv <- if (length(ph) >= 2) rbindlist(lapply(names(ph[[length(ph)]]), function(nm)
  data.table(block = nm, last_change = mean(abs(ph[[length(ph)]][[nm]] - ph[[length(ph) - 1]][[nm]]))))) else data.table(block = NA_character_, last_change = NA_real_)
fwrite(conv, file.path(FIT_DIR, "convergence_stageB.csv"))

writeLines(c(
  sprintf("# FINN-model-FIA (oregon), zoo v0.2 recipe, stage B (%s)", PARAMS$prep_version), "",
  sprintf("Fitted %s by `calibration/v02/recipe/oregon/fit_stageB.R` in %.1f min (%.1f s per epoch) on %s with %d thread(s): FINN %s from `%s`, torch %s, optimizer optim_adam, lr %g with %s decay.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), fit_minutes, 60 * fit_minutes / PARAMS$epochs, PARAMS$device, PARAMS$threads,
          as.character(packageVersion("FINN")), basename(dirname(find.package("FINN"))), as.character(packageVersion("torch")), PARAMS$lr, PARAMS$lr_scheduler), "",
  sprintf("Data: `prepared/%s/%s` (%d sites, %d species, %d steps of %d yr). Holdout for validation: `%s`. Frozen demography: `model_stageA.pt` (competition, growth, mortality incl. the U-shape and both networks, dropout %g while frozen); max |change| of the frozen blocks %s.",
          PARAMS$prep_version, PARAMS$train_set, uniqueN(obs$siteID), n_species, nrow(iv), PARAMS$step_years, PARAMS$holdout_set, PARAMS$frozen_dropout,
          paste(sprintf("%s %.1g", names(dev), dev), collapse = ", ")), "",
  sprintf("Regeneration: FINN::regeneration_saturation (light x exp(linear env on `%s`), Beverton-Holt cap K per species, start %g/ha/step, reg_floor %g, recruits drawn with the \"%s\" sampler), intercept warm-started at the observed mean recruitment per step; theta per species = intercept-only ML on the observed recruit counts (reg x %.5f ha), clamped to [%g, %g], frozen: %s.",
          paste(deparse(env_formula), collapse = ""), PARAMS$reg_K_init, PARAMS$reg_floor, PARAMS$reg_sampler, w_obs, PARAMS$theta_min, PARAMS$theta_max,
          paste(sprintf("%s %.3g", REG$species_name, REG$theta), collapse = ", ")), "",
  sprintf("Fit: %d epochs, batch %d, weights %s (weight_scale %s), losses %s, growth and mortality at the site x species level (no tree table; period_mode per_site), seed %d.",
          PARAMS$epochs, PARAMS$batchsize, PARAMS$weights, PARAMS$weight_scale,
          paste(sprintf("%s=%s", names(PARAMS$loss_family), PARAMS$loss_family), collapse = " "), PARAMS$seed), "",
  sprintf("Baseline-scaled losses (1 = no better than the intercept), first and last epoch: %s.",
          paste(sprintf("%s %.3f -> %.3f", names(loss_hist)[1:6], as.numeric(loss_hist[1, 1:6]), as.numeric(loss_hist[.N, 1:6])), collapse = ", ")), "",
  sprintf("Predicted / observed recruitment on the training plots (fit windows): %.2f at the start, %.2f after stage B.", ratio_start, ratio_end), "",
  sprintf("Fitted cap K (recruits/ha/step; observed max per-ha rate per step of the species in brackets): %s (`regeneration_parameters.csv`). K ran away if it exceeds the observed maximum by far (v0.1: up to 131,000); the pilots' alternative is a fixed cap at 2 x the 95th percentile.",
          paste(sprintf("%s %.1f (%.0f)", REG$species_name, REG$K_fit, REG$obs_max), collapse = ", ")), "",
  sprintf("Largest parameter changes between the last two snapshots (every %d epochs): %s.", PARAMS$checkpoints,
          paste(sprintf("%s %.2g", conv[order(-last_change)][1:5, block], conv[order(-last_change)][1:5, last_change]), collapse = ", ")), "",
  "Validation report: `Rscript calibration/v02/recipe/oregon/validation.R`."
), file.path(FIT_DIR, "fit_summary_stageB.md"))
message(sprintf("stage B done in %.1f min (%.1f s/epoch); recruitment ratio %.2f -> %.2f; K range %.1f-%.1f; final losses: %s", fit_minutes, 60 * fit_minutes / PARAMS$epochs,
                ratio_start, ratio_end, min(REG$K_fit), max(REG$K_fit), paste(sprintf("%s=%.3f", names(loss_hist)[1:6], as.numeric(loss_hist[.N, 1:6])), collapse = " ")))
