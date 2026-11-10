#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/pilot_oregon_O0_v01.R
##
## WHAT   Pilot O0, the REFERENCE: the v0.1 Oregon model settings (calibration/
##        oregon/02_fit.R) on the v0.2 DATA (oregon-10yr-v2: 18 species under the
##        agreed rare-species rule, no `other`), with the pilot budget (300 epochs
##        instead of 500, one thread, same seed), scored with the pilot
##        evaluation (holdout against persistence, change metrics, mortality
##        AUC, recruitment ratio, 400-year bare-ground run, errors stratified by
##        the end-census out-of-S share). Every other pilot differs from this
##        one by the settings named in its header; the data are the same.
## DATA   oregon-10yr-v2/train_blocked (fit), holdout_blocked (evaluation)
## RUN    Rscript calibration/v02/pilot_oregon_O0_v01.R   (repo root, ~15 min, one thread)
## OUTPUT calibration/v02/pilots/oregon_O0_v01/{metrics.csv, summary.csv, predictions.csv, bare_ground.csv, loss_history.csv}
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
source(file.path("calibration", "v02", "pilot_eval.R"))
torch_set_num_threads(1L)
## torch >= 0.14 has the faster ignite Adam (FINN's default); the cluster's torch 0.13 falls back to optim_adam
OPTIMIZER <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
DATA_DIR <- normalizePath(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))[
  dir.exists(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia")))][1])
LABEL <- "oregon_O0_v01"; OUT <- file.path("calibration", "v02", "pilots", LABEL); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  prep_version = "oregon-10yr-v2", train_set = "train_blocked", holdout_set = "holdout_blocked",   # v0.2 data: S = 18 species, < 2 % out-of-S site rule (every pilot uses it)
  patch_size_ha = 0.01681, n_patches = 4L, recruits_dbh = 12.9,
  env_formula = ~ temp + prec, hidden = NULL,
  reg_K_init = 500, reg_floor = 0.01,            # v0.1: K trainable, init 500 per ha per decade
  epochs = 300L, batchsize = 50L, lr = 0.01, lr_scheduler = "none",   # v0.1 settings at the pilot budget
  weights = "auto", weight_scale = "auto",
  loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  bare_steps = 40L, bare_patches = 10L, seed = 42L
)

train <- pe_read_set(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set))
hold  <- pe_read_set(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$holdout_set))
n_species <- nrow(train$species)

FINN.seed(PARAMS$seed)
model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(PARAMS$env_formula, FINN::growth, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(PARAMS$env_formula, FINN::mortality, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(PARAMS$env_formula, FINN::regeneration_saturation, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       custom_parameters = list(reg_logK = rep(log(PARAMS$reg_K_init), n_species))))
t0 <- Sys.time()
fit(model, env = train$env, data = train$obs, init_cohort = makeInitCohorts(train$init, Nspecies = n_species),
    epochs = PARAMS$epochs, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler,
    weights = PARAMS$weights, weight_scale = PARAMS$weight_scale, loss_family = PARAMS$loss_family,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, env_autoscale = TRUE, device = "cpu",
    optimizer = OPTIMIZER, folder = NULL, checkpoints = 10L, plot_progress = FALSE)
minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
lh <- as.data.table(do.call(rbind, model$history)); setnames(lh, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(lh))]); lh[, epoch := .I]
fwrite(lh, file.path(OUT, "loss_history.csv"))
torch_save(lapply(model$parameters, function(p) p$detach()$cpu()$clone()), file.path(DATA_DIR, "fits", paste0("pilot_", LABEL, ".pt")))

## bare ground: three climate points of the training set, 40 steps = 400 years
se <- train$env[, .(temp = mean(temp), prec = mean(prec)), by = siteID][order(temp)]
bare <- se[round(quantile(seq_len(.N), c(0.1, 0.5, 0.9))), .(siteID = .I, temp, prec)]
res <- pilot_evaluate(model, hold, LABEL, OUT, n_patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, step_years = 10L,
                      bare_env = bare, bare_years = PARAMS$bare_steps, bare_patches = PARAMS$bare_patches)
pilot_stratify(res$predictions, hold$sites, OUT)
K <- exp(as.numeric(as_array(model$reg_logK))); theta <- as.numeric(as_array(model$par_theta_recruits))
res$summary[, `:=`(minutes = round(minutes, 1), epochs = PARAMS$epochs, K_median = median(K), K_max = max(K), theta_median = median(theta),
                   loss_last = paste(sprintf("%s %.3f", names(lh)[1:6], as.numeric(lh[.N, 1:6])), collapse = " "))]
fwrite(res$summary, file.path(OUT, "summary.csv"))
message(sprintf("%s done in %.0f min; K median %.1f max %.0f; theta median %.3f", LABEL, minutes, median(K), max(K), median(theta)))
