#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-zoo | FINN-model-FIA/calibration/v02/recipe/conus/fit_stageA.R
##
## WHAT   Stage A of the zoo v0.2 recipe for the conterminous-US variant (annual
##        step, per-plot intervals of 4-12 years, 36 species, data conus-annual-v2):
##        competition, growth and mortality are fitted with the OBSERVED
##        recruitment fed into the simulation, so that no regeneration parameter
##        is trained here and the demography is not asked to explain recruits.
##        Stage B (fit_stageB.R) then fits the regeneration process alone.
##        Copied from the reference implementation
##        FINN-model-IFN-Spain/calibration/v02/recipe/fit_stageA.R (cb41778); the
##        per-model block of PARAMS and the data handling differ as stated below.
## WHY    The recipe (FINN-zoo-data/RECIPE_v02.md, decided 2026-11-10):
##        - mortality: FINN::mortality's linear predictor plus a U-shaped size term
##          sen * (dbh/100)^2 (custom function `mortality_ushape`, no FINN change);
##        - environment: a SELU network (hidden 20-20, dropout 0.1) on temperature,
##          precipitation and the water deficit for growth and for mortality; no
##          coordinates, so the tool transfers to new sites;
##        - recruitment: the observed recruits of every plot x species x interval,
##          fed as a per-patch rate over the years of the interval (custom function
##          `regeneration_observed`, env columns reg_obs_<s>), frozen;
##        - losses: all state and rate terms (dbh, ba, trees, growth, mortality);
##          growth and mortality scored per tracked tree; no regeneration term
##          (it is observed here); weights "auto";
##        - cosine learning-rate decay, plain Adam, one thread, light_mask TRUE,
##          FINN's default mortality floor and gamma draw.
##        DEVIATIONS from the reference, with their reason (data facts of this variant):
##        1. The observed recruits are NOT placed as whole trees per patch and year.
##           A plot is 4 patches of 0.01681 ha with one site-level env, so one
##           observed recruit is 0.25 trees per patch; FINN rounds every patch and
##           year to whole trees, which drops 1-recruit cells (54 % of the recruiting
##           cells of the training set) and keeps 67 % of the recruits in total. The
##           schedule therefore carries the per-patch MEAN rate (reg / p per ha per
##           year on the rows of the interval) and FINN draws the recruits with
##           `regeneration_sampler = "refined"` (exact torch_poisson below a mean of 2
##           per patch and year; the default "gaussian" sampler delivers ~0 % of the
##           means this data set has, 0.0025-0.1 per patch and year) at a near-Poisson
##           dispersion (theta 100, frozen). The fed recruitment equals the observed
##           per cell in expectation; the check below verifies the scheduled mean
##           (ratio 1.0000) and the sampler's mean at the smallest rate of the data.
##           The same sampler is used in stage B and the simulation model.
##        2. Several censuses per plot (1-4; year != period_length): the schedule and
##           the fed-recruitment check work per interval (env rows y_end - p + 1 .. y_end).
##        3. The climate columns are the YEARLY TerraClimate values of
##           conus-annual-v2 (tavg_c, prec_mm, cwd_z of every simulated year), not
##           static normals: the brief names this data version and its water deficit
##           for the conus variant, and the version carries no normals.
## DATA   ../../FINN-zoo-data/fia/prepared/conus-annual-v2/train
##        (obs_dt, env_dt, init_cohort_dt, tree_obs_dt, species_dt)
## NEEDS  FINN 0.2.0 (FINNverse/FINN main, commit 0848b60) and torch; the library
##        ~/R/lib-finn-0848b60 is put first on .libPaths() when it exists (laptop and
##        cluster), otherwise R_LIBS_USER must point at it.
## RUN    Rscript calibration/v02/recipe/conus/fit_stageA.R      (from the repo root)
##        cluster: hpc/submit_recipe_conus.sh (stage A, then stage B and validation chained)
## OUTPUT <FIT_DIR>/model_stageA.pt          the training model (observed recruitment fed)
##        <FIT_DIR>/loss_history_stageA.csv  one row per epoch
##        <FIT_DIR>/species_parameters_stageA.csv, mortality_size_curve_stageA.csv, convergence_stageA.csv
##        <FIT_DIR>/fit_summary_stageA.md
##        FIT_DIR = ../../FINN-zoo-data/fia/fits/conus-annual-v2-recipe (outside git)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
if (dir.exists("~/R/lib-finn-0848b60")) .libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({
  library(FINN)
  library(data.table)
  library(torch)
})
stopifnot(packageVersion("FINN") == "0.2.0")

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "fia")
FIT_DIR  <- file.path(DATA_DIR, "fits", "conus-annual-v2-recipe")
dir.create(FIT_DIR, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  ## ---- per model: what changes from data set to data set ----
  prep_version  = "conus-annual-v2",       # v0.2 data: 36 species, no "other", start-census site rule x = 2 %, complete ingrowth target
  train_set     = "train",                 # 1,500 training plots in the retained 2-degree blocks
  env_vars      = c("tavg_c", "prec_mm", "cwd_z"),  # yearly TerraClimate: mean temperature (C), precipitation (mm), growing-season water-deficit z-anomaly; the network learns the interactions
  patch_size_ha = 0.01681,                 # one FINN patch = one FIA subplot (7.3 m radius)
  n_patches     = 4L,                      # the four subplots of a plot
  recruits_dbh  = 12.9,                    # recruits enter just above the 12.7 cm (5 in) tally threshold
  reg_sampler   = "refined",               # deviation 1: exact Poisson draw of the per-patch recruits below a mean of 2 (FINN's default "gaussian" loses the small per-patch means of 0.0168-ha patches)
  reg_obs_theta = 100,                     # near-Poisson dispersion of the fed recruitment draw (frozen; the observed rate is the mean)
  ## ---- the recipe (same in every model) ----
  hidden        = c(20L, 20L),             # environment network of growth and mortality (SELU, as the zoo fits)
  dropout       = 0.1,                     # as the zoo fits; active while training, off in simulateForest (eval mode)
  sen_max       = 5,                       # U-shape cap: at most +5 on the mortality logit at dbh 1 m (BWI final model)
  sen_init      = 0.05,                    # U-shape start: +0.05 on the logit at 1 m, practically a linear size term
  growth_scale  = 0.01,                    # growth enters mortality as growth / 0.01 (growth in percent; FINN 0.2.0's default)
  reg_floor     = 1e-3,                    # FINN's additive recruitment floor (stems/ha/yr); unused by the observed recruitment, kept for stage B
  reg_obs_eps   = 1e-8,                    # recruits/ha/yr standing in for a year without observed recruits (log scale)
  light_mask    = TRUE,                    # 0.2.0's masked competition: padding slots do not set the height-quantile edges
  n_quantiles   = 10L,                     # FINN's default height-quantile competition
  epochs        = 300L,                    # cosine decay to lr/100 over the epochs; the pilots settled within 200
  batchsize     = 100L,                    # as the zoo pilots
  lr            = 0.005,                   # as the zoo pilots (0.01 ran away in the mortality term on FIA)
  lr_scheduler  = "cosine",                # the recipe: settles reproducibly
  weights       = "auto",                  # each loss scaled by its intercept-only baseline (FINN 0.2.0); no numeric vector needed, see loss_family
  weight_scale  = "auto",                  # stated for the card: no numeric weights, so the scale setting is moot
  loss_family   = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse",
                    mortality = "binomial", regeneration = "none"),   # all state and rate terms; recruitment is observed here, so no term
  loss_aggregation = list(growth = "tree", mortality = "tree"),   # per tracked tree (tree_obs_dt); needs period_mode per_site
  checkpoints   = 10L,                     # a parameter snapshot every 10 epochs (folder NULL: no model files)
  threads       = 1L,                      # one thread: several zoo fits share a machine
  device        = "cpu",
  seed          = 1L
)
torch_set_num_threads(PARAMS$threads)
stopifnot(identical(formals(FINN::mortality)$growth_scale, PARAMS$growth_scale))   # the U-shape kernel below uses FINN's growth scale
PS <- PARAMS$patch_size_ha

## ---- 1. data -----------------------------------------------------------------
train_dir <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)
obs   <- fread(file.path(train_dir, "obs_dt.csv"))         # site x species at every census year (interval of period_length years)
env   <- fread(file.path(train_dir, "env_dt.csv"))[, c("siteID", "year", PARAMS$env_vars), with = FALSE]
init  <- fread(file.path(train_dir, "init_cohort_dt.csv")) # the first-census stand
trees <- fread(file.path(train_dir, "tree_obs_dt.csv"))    # one row per tracked tree and interval (growth, died)
species_dt <- fread(file.path(train_dir, "species_dt.csv"))
n_species <- nrow(species_dt)
stopifnot(max(obs$species) == n_species, all(obs$year >= obs$period_length), uniqueN(init$patchID) == PARAMS$n_patches)
init_cohort <- makeInitCohorts(init, Nspecies = n_species)
Tmax <- max(env$year)
## the intervals of a plot follow each other without overlap (the schedule below sums per row)
iv <- unique(obs[, .(siteID, y_end = year, y0 = year - period_length)])[order(siteID, y_end)]
stopifnot(iv[, all(y0 == shift(y_end, fill = 0L)), by = siteID][, all(V1)], max(iv$y_end) <= Tmax)
message(sprintf("train: %d sites, %d species, %d intervals of %d..%d yr (1-%d per site), Tmax %d, %d tracked tree intervals (%d died)",
                uniqueN(obs$siteID), n_species, nrow(iv), min(obs$period_length), max(obs$period_length), iv[, .N, by = siteID][, max(N)],
                Tmax, nrow(trees), sum(trees$died)))

## The observed recruitment as env columns reg_obs_<s>. The observed recruits of a
## plot x species x interval, n = reg x (4 x patch area) trees on the plot, are a
## per-patch mean of n / 4 over the p years of the interval; the rows y_end - p + 1
## .. y_end carry the per-ha rate reg / p, the rows outside every interval none.
## The column holds log(rate + eps), which `regeneration_observed` turns back into
## the per-ha rate FINN expects; FINN draws whole recruits per patch and year from
## it with the "refined" sampler (deviation 1 in the header).
reg_cols <- paste0("reg_obs_", seq_len(n_species))
sched <- obs[, .(siteID, species, y_end = year, p = period_length, n = round(pmax(reg, 0) * PARAMS$n_patches * PS), rate = pmax(reg, 0) / period_length)]
stopifnot(max(abs(obs$reg * PARAMS$n_patches * PS - sched$n)) < 0.01)   # whole recruits on the plot
sched <- sched[, .(year = seq_len(Tmax)), by = .(siteID, species, y_end, p, n, rate)]
sched[, inside := year > y_end - p & year <= y_end]
sched[, fed_rate := fifelse(inside, rate, 0)]
REG <- sched[, .(fed_rate = sum(fed_rate)), by = .(siteID, species, year)]   # the intervals do not overlap (checked above)
REG[, col := paste0("reg_obs_", species)]
REG <- dcast(REG, siteID + year ~ col, value.var = "fed_rate")
setcolorder(REG, c("siteID", "year", reg_cols))
for (cc in reg_cols) set(REG, j = cc, value = log(REG[[cc]] + PARAMS$reg_obs_eps))
env_obsreg <- merge(env, REG, by = c("siteID", "year"))
setorder(env_obsreg, siteID, year)
stopifnot(nrow(env_obsreg) == nrow(env))
n_fed <- sched[year == 1, sum(n)]
message(sprintf("observed recruits fed: %d trees on %d plots (%.1f %% of the site x species x interval cells recruit)",
                n_fed, sched[year == 1 & n > 0, uniqueN(siteID)], 100 * mean(obs$reg > 0)))
## the sampler must deliver the scheduled mean at the smallest per-patch rate of the
## data (one recruit per plot over the longest interval): the default "gaussian"
## sampler would not (it rounds lambda + eps * sqrt(lambda) and delivers ~0 % here)
lam_min <- 1 / PARAMS$n_patches / max(obs$period_length)
FINN.seed(PARAMS$seed)
draw <- as.numeric(FINN:::rnbinom_torch_refined(torch_full(c(200000L), lam_min), torch_tensor(PARAMS$reg_obs_theta)))
sampler_ratio <- mean(round(draw)) / lam_min
stopifnot(abs(sampler_ratio - 1) < 0.05)
message(sprintf("refined sampler at the smallest per-patch mean %.4f: delivered / expected %.3f", lam_min, sampler_ratio))

## ---- 2. the custom process functions (the recipe's code blocks) --------------------
## U-shaped mortality: FINN::mortality's linear predictor (light, size, growth in
## percent) plus sen * (dbh/100)^2 with sen = sen_max * sigmoid(mort_sen_raw) per
## species in [0, sen_max]. Custom parameters: mort_sen_raw (trainable), mort_sen_max
## (fixed). The growth_scale formal is the one FINN::mortality has, so FINN
## initialises the growth coefficient on the same scale.
mortality_ushape <- function(dbh, species, trees, parMort, pred, light, base_steepness = 5, debug = FALSE, growth = NULL, growth_scale = 0.01) {
  if (is.null(growth)) growth <- self$g
  d   <- dbh / 100
  sen <- (self$mort_sen_max * torch::torch_sigmoid(self$mort_sen_raw))[species]
  torch::torch_sigmoid(pred + growth / growth_scale * parMort[, 3][species] + light * parMort[, 1][species] +
                       parMort[, 2][species] * d + sen * d * d)
}
sen_raw_init <- rep(qlogis(PARAMS$sen_init / PARAMS$sen_max), n_species)
## Observed recruitment: the env design matrix ~ 0 + reg_obs_1 + ... + reg_obs_S passes
## through an identity linear layer (initEnv = diag, frozen), so pred[, s] is the
## log per-ha rate of species s in that year; the function returns exp(pred) for
## every patch and ignores light and the species parameters.
regeneration_observed <- function(species, parReg, pred, light) {
  torch_exp(pred)[, NULL]$`repeat`(c(1, species$shape[2], 1))
}
obsreg_formula <- as.formula(paste("~ 0 +", paste(reg_cols, collapse = " + ")))
env_formula    <- as.formula(paste("~", paste(PARAMS$env_vars, collapse = " + ")))
## FINN z-scales the env columns it is given (env_autoscale); the reg_obs columns
## must stay as they are, so the scaling is named per column: climate "auto",
## recruitment "identity".
autoscale_spec <- c(setNames(rep("auto", length(PARAMS$env_vars)), PARAMS$env_vars), setNames(rep("identity", n_species), reg_cols))

## ---- 3. the training model, written out ---------------------------------------
FINN.seed(PARAMS$seed)
model <- finn(
  N_species    = n_species,
  recruits_dbh = PARAMS$recruits_dbh,
  reg_floor    = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE,
                                       n_quantiles = PARAMS$n_quantiles, light_mask = PARAMS$light_mask),
  growth_process       = createProcess(env_formula, FINN::growth,
                                       hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(env_formula, mortality_ushape,
                                       hidden = PARAMS$hidden, dropout = PARAMS$dropout,
                                       optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       custom_parameters = list(mort_sen_raw = list(init = sen_raw_init, trainable = TRUE),
                                                                mort_sen_max = list(init = PARAMS$sen_max, trainable = FALSE))),
  regeneration_process = createProcess(obsreg_formula, regeneration_observed,
                                       initEnv = list(diag(n_species)),
                                       dispersion_parameter = rep(PARAMS$reg_obs_theta, n_species),
                                       optimizeSpecies = FALSE, optimizeEnv = FALSE,
                                       sample_regeneration = TRUE, regeneration_sampler = PARAMS$reg_sampler)
)
## the dispersion of the fed draw is not a parameter of this stage (the draw is
## reparameterised, so a trainable theta would receive state-loss gradients)
invisible(model$par_theta_recruits_raw$requires_grad_(FALSE))
stopifnot(identical(model$regeneration_sampler, PARAMS$reg_sampler), isTRUE(model$sample_regeneration))
## FINN scores the per-tree losses on its per-site path only (chosen automatically
## when the intervals differ between sites; forced here so a data set with one
## interval length is scored the same way)
model$period_mode <- "per_site"

## ---- 4. fit ----------------------------------------------------------------------
FINN.seed(PARAMS$seed)
t0 <- Sys.time()
fit(model,
    env           = env_obsreg,
    data          = obs,
    tree_data     = trees,
    init_cohort   = init_cohort,
    epochs        = PARAMS$epochs,
    batchsize     = PARAMS$batchsize,
    lr            = PARAMS$lr,
    lr_scheduler  = PARAMS$lr_scheduler,
    optimizer     = torch::optim_adam,     # plain Adam, the same on the laptop and the cluster
    weights       = PARAMS$weights,
    weight_scale  = PARAMS$weight_scale,
    loss_family   = PARAMS$loss_family,
    loss_aggregation = PARAMS$loss_aggregation,
    patches       = PARAMS$n_patches,
    patch_size    = PS,
    env_autoscale = autoscale_spec,
    device        = PARAMS$device,
    folder        = NULL,
    checkpoints   = PARAMS$checkpoints,
    plot_progress = FALSE)
fit_minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

## save first, summarise later
torch_save(model, file.path(FIT_DIR, "model_stageA.pt"))
loss_hist <- as.data.table(do.call(rbind, model$history))
setnames(loss_hist, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(loss_hist))])
loss_hist[, epoch := .I]
fwrite(loss_hist, file.path(FIT_DIR, "loss_history_stageA.csv"))

## ---- 5. checks and what was fitted ---------------------------------------------------
## the fed recruitment (the scheduled mean, r_mean_ha) must equal the observed
## recruits over fit()'s windows: the observation at year y_end is driven by env
## rows y_end - p + 1 .. y_end, and simulateForest() labels the year driven by env
## row k as year k (checked: the rate scheduled on row 1 appears in output year 1)
model$eval()
FINN.seed(PARAMS$seed)
S <- as.data.table(simulateForest(model, env = env_obsreg, init_cohort = init_cohort, patches = PARAMS$n_patches,
                                  patch_size = PS, device = "cpu")$wide$site)
for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
W <- iv[, .(year = seq.int(y0 + 1L, y_end)), by = .(siteID, y_end)]
fed <- merge(W, S[, .(siteID, year, species, r_mean_ha)], by = c("siteID", "year"), allow.cartesian = TRUE)[
  , .(fed = sum(r_mean_ha) * PARAMS$n_patches * PS), by = .(siteID, y_end, species)]
fed <- merge(fed, sched[year == 1, .(siteID, y_end, species, n)], by = c("siteID", "y_end", "species"))
fed_ratio <- sum(fed$fed) / sum(fed$n)
stopifnot(abs(fed_ratio - 1) < 1e-3, max(abs(fed$fed - fed$n)) < 0.02)
## the fitted species parameters and the U-shape per species
spar <- function(block, nm) `colnames<-`(as.data.table(as.matrix(model[[block]]$detach()$cpu())), nm)
SPP <- cbind(species_dt[order(species)], spar("par_competition", c("comp_height", "comp_strength")),
             spar("par_growth", c("growth_light", "growth_size")), spar("par_mortality", c("mort_light", "mort_size", "mort_growth")))
fwrite(SPP, file.path(FIT_DIR, "species_parameters_stageA.csv"))
sen_raw <- as.numeric(model$mort_sen_raw$detach()$cpu())
SEN <- data.table(species = seq_len(n_species), p_dbh = SPP$mort_size, sen = PARAMS$sen_max * plogis(sen_raw))
SEN[, dbh_min_mortality_cm := fifelse(sen > 1e-3 & p_dbh < 0, -p_dbh / (2 * sen) * 100, Inf)]   # dbh of minimum mortality
SEN <- merge(species_dt, SEN, by = "species")
fwrite(SEN, file.path(FIT_DIR, "mortality_size_curve_stageA.csv"))
## convergence: mean |change| of every block between the last two snapshots
ph <- model$param_history
conv <- if (length(ph) >= 2) rbindlist(lapply(names(ph[[length(ph)]]), function(nm)
  data.table(block = nm, last_change = mean(abs(ph[[length(ph)]][[nm]] - ph[[length(ph) - 1]][[nm]]))))) else data.table(block = NA_character_, last_change = NA_real_)
fwrite(conv, file.path(FIT_DIR, "convergence_stageA.csv"))

writeLines(c(
  sprintf("# FINN-model-FIA (conus), zoo v0.2 recipe, stage A (%s)", PARAMS$prep_version), "",
  sprintf("Fitted %s by `calibration/v02/recipe/conus/fit_stageA.R` in %.1f min (%.1f s per epoch) on %s with %d thread(s): FINN %s from `%s`, torch %s, optimizer optim_adam, lr %g with %s decay.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), fit_minutes, 60 * fit_minutes / PARAMS$epochs, PARAMS$device, PARAMS$threads,
          as.character(packageVersion("FINN")), basename(dirname(find.package("FINN"))), as.character(packageVersion("torch")), PARAMS$lr, PARAMS$lr_scheduler), "",
  sprintf("Data: `prepared/%s/%s` (%d sites, %d species, %d intervals of %d..%d yr, %d tracked tree intervals of which %d died).",
          PARAMS$prep_version, PARAMS$train_set, uniqueN(obs$siteID), n_species, nrow(iv), min(obs$period_length), max(obs$period_length), nrow(trees), sum(trees$died)), "",
  sprintf("Model: competition (light_mask %s, %d quantiles); growth FINN::growth and mortality `mortality_ushape` (sen in [0, %g], start %g at 1 m; growth / %g), each with a SELU network hidden = c(%s), dropout %g on `%s` (yearly values); regeneration = the observed recruits fed as a per-patch mean rate per interval (`regeneration_observed`, %d trees on %d plots; fed / observed over the fit windows %.4f) and drawn with FINN's \"%s\" sampler at theta %g (delivered / expected at the smallest per-patch mean %.4f: %.3f).",
          PARAMS$light_mask, PARAMS$n_quantiles, PARAMS$sen_max, PARAMS$sen_init, PARAMS$growth_scale, paste(PARAMS$hidden, collapse = ", "), PARAMS$dropout,
          paste(deparse(env_formula), collapse = ""), n_fed, sched[year == 1 & n > 0, uniqueN(siteID)], fed_ratio, PARAMS$reg_sampler, PARAMS$reg_obs_theta, lam_min, sampler_ratio), "",
  sprintf("Fit: %d epochs, batch %d, weights %s (weight_scale %s), losses %s, growth and mortality per tracked tree (period_mode per_site), env_autoscale climate auto / reg_obs identity, %d patches of %.5f ha, recruits_dbh %g, seed %d.",
          PARAMS$epochs, PARAMS$batchsize, PARAMS$weights, PARAMS$weight_scale,
          paste(sprintf("%s=%s", names(PARAMS$loss_family), PARAMS$loss_family), collapse = " "),
          PARAMS$n_patches, PS, PARAMS$recruits_dbh, PARAMS$seed), "",
  sprintf("Baseline-scaled losses (1 = no better than the intercept), first and last epoch: %s.",
          paste(sprintf("%s %.3f -> %.3f", names(loss_hist)[1:5], as.numeric(loss_hist[1, 1:5]), as.numeric(loss_hist[.N, 1:5])), collapse = ", ")), "",
  sprintf("U-shape: sen %s; dbh of minimum mortality (cm): %s (`mortality_size_curve_stageA.csv`).",
          paste(sprintf("%.3f", range(SEN$sen)), collapse = "-"),
          paste(sprintf("%s %.0f", SEN$species_name, SEN$dbh_min_mortality_cm), collapse = ", ")), "",
  sprintf("Largest parameter changes between the last two snapshots (every %d epochs): %s.", PARAMS$checkpoints,
          paste(sprintf("%s %.2g", conv[order(-last_change)][1:5, block], conv[order(-last_change)][1:5, last_change]), collapse = ", ")), "",
  "Deviations from the reference implementation: (1) recruits fed as a per-patch mean rate and drawn with the refined sampler (4 x 0.0168 ha patches, see the script header); (2) per-interval windows for plots with several censuses; (3) yearly TerraClimate columns instead of static normals (the data version carries no normals).", "",
  "Next: `Rscript calibration/v02/recipe/conus/fit_stageB.R` (regeneration), then `validation.R`."
), file.path(FIT_DIR, "fit_summary_stageA.md"))
message(sprintf("stage A done in %.1f min (%.1f s/epoch); final losses: %s", fit_minutes, 60 * fit_minutes / PARAMS$epochs,
                paste(sprintf("%s=%.3f", names(loss_hist)[1:5], as.numeric(loss_hist[.N, 1:5])), collapse = " ")))
