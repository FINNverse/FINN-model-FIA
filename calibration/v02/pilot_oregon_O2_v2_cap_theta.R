#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/pilot_oregon_O2_v2_cap_theta.R
##
## WHAT   Pilot O2 = O1 (fixed data-derived recruitment cap, cosine lr) with
##        two more changes, each tied to an audit finding:
##        1. DATA v2 (oregon-10yr-v2, calibration/v02/prepare_oregon_v2.R): the
##           lumped `other` class is dropped (audit G: 1 % of the observed
##           basal area and recruits, 27 % of the v0.1 predicted recruits).
##        2. FROZEN RECRUITMENT DISPERSION (brief trap 7): the per-species NB
##           theta is set from the v0.1 production fit (its ML-like values,
##           median 0.08) clipped to [0.05, 2] and not trained, so the
##           regeneration loss cannot flatten itself by letting theta fall.
##        Everything else as O1.
## RUN    Rscript calibration/v02/pilot_oregon_O2_v2_cap_theta.R   (repo root, ~15 min, one thread)
## OUTPUT calibration/v02/pilots/oregon_O2_v2_cap_theta/
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
source(file.path("calibration", "v02", "pilot_eval.R"))
source(file.path("R", "definition.R"))
torch_set_num_threads(1L)
DATA_DIR <- normalizePath(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))[
  dir.exists(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia")))][1])
LABEL <- "oregon_O2_v2_cap_theta"; OUT <- file.path("calibration", "v02", "pilots", LABEL); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  prep_version = "oregon-10yr-v2", train_set = "train_blocked", holdout_set = "holdout_blocked",
  patch_size_ha = 0.01681, n_patches = 4L, recruits_dbh = 12.9,
  env_formula = ~ temp + prec, hidden = NULL,
  cap_quantile = 0.95, cap_factor = 2, cap_floor = 20,
  reg_floor = 0.01,
  theta_source = "oregon",                 # the v0.1 production fit's per-species theta (its first 27 species)
  theta_bounds = c(0.05, 2),               # clipped, then frozen
  epochs = 300L, batchsize = 50L, lr = 0.01, lr_scheduler = "cosine",
  weights = "auto", weight_scale = "auto",
  loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  bare_steps = 40L, bare_patches = 10L, seed = 42L
)

train <- pe_read_set(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set))
hold  <- pe_read_set(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$holdout_set))
n_species <- nrow(train$species)
K_cap <- pe_reg_cap(train$obs, train$species, q = PARAMS$cap_quantile, factor = PARAMS$cap_factor, floor_cap = PARAMS$cap_floor, per_step = 1)
## theta of the v0.1 fit (28 species incl. other; v2 keeps species 1..27 with the same codes)
m01 <- build_model(PARAMS$theta_source, torch_load(file.path(DATA_DIR, "fits", PARAMS$theta_source, "parameters.pt")),
                   env_scaling = as.data.frame(fread(file.path("calibration", PARAMS$theta_source, "env_scaling.csv"))))
theta0 <- pmin(pmax(as.numeric(as_array(m01$par_theta_recruits))[seq_len(n_species)], PARAMS$theta_bounds[1]), PARAMS$theta_bounds[2])
message("theta (frozen): ", paste(sprintf("%s %.3f", train$species$species_name, theta0), collapse = ", "))

FINN.seed(PARAMS$seed)
model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(PARAMS$env_formula, FINN::growth, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(PARAMS$env_formula, FINN::mortality, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(PARAMS$env_formula, FINN::regeneration_saturation, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       dispersion_parameter = theta0,
                                       custom_parameters = list(reg_logK = list(init = log(K_cap), trainable = FALSE))))
model$par_theta_recruits_raw$requires_grad_(FALSE)      # frozen dispersion (brief trap 7)
t0 <- Sys.time()
fit(model, env = train$env, data = train$obs, init_cohort = makeInitCohorts(train$init, Nspecies = n_species),
    epochs = PARAMS$epochs, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler,
    lr_scheduler_params = list(T_max = PARAMS$epochs),
    weights = PARAMS$weights, weight_scale = PARAMS$weight_scale, loss_family = PARAMS$loss_family,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, env_autoscale = TRUE, device = "cpu",
    folder = NULL, checkpoints = 10L, plot_progress = FALSE)
minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
lh <- as.data.table(do.call(rbind, model$history)); setnames(lh, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(lh))]); lh[, epoch := .I]
fwrite(lh, file.path(OUT, "loss_history.csv"))
torch_save(lapply(model$parameters, function(p) p$detach()$cpu()$clone()), file.path(DATA_DIR, "fits", paste0("pilot_", LABEL, ".pt")))

se <- train$env[, .(temp = mean(temp), prec = mean(prec)), by = siteID][order(temp)]
bare <- se[round(quantile(seq_len(.N), c(0.1, 0.5, 0.9))), .(siteID = .I, temp, prec)]
res <- pilot_evaluate(model, hold, LABEL, OUT, n_patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, step_years = 10L,
                      bare_env = bare, bare_years = PARAMS$bare_steps, bare_patches = PARAMS$bare_patches)
theta <- as.numeric(as_array(model$par_theta_recruits))
res$summary[, `:=`(minutes = round(minutes, 1), epochs = PARAMS$epochs, K_median = median(K_cap), K_max = max(K_cap), theta_median = median(theta),
                   loss_last = paste(sprintf("%s %.3f", names(lh)[1:6], as.numeric(lh[.N, 1:6])), collapse = " "))]
fwrite(res$summary, file.path(OUT, "summary.csv"))
message(sprintf("%s done in %.0f min; theta median %.3f (frozen)", LABEL, minutes, median(theta)))
