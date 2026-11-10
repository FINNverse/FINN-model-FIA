#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/pilot_conus_C0_v01.R
##
## WHAT   Pilot C0, the REFERENCE: the v0.1 CONUS model (calibration/conus/
##        02_fit.R) on the fixed PILOT SUBSET (300 training plots drawn from
##        conus-annual-v1/train with seed 7, 200 holdout plots from
##        conus-annual-v1/holdout with seed 8), 250 epochs, one thread, scored
##        with the pilot evaluation (holdout against persistence, change
##        metrics, per-tree AUC and growth R2, recruitment ratio, 300-year
##        bare-ground run). The other conus pilots use the same plots.
## RUN    Rscript calibration/v02/pilot_conus_C0_v01.R   (repo root, ~40 min, one thread)
## OUTPUT calibration/v02/pilots/conus_C0_v01/
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
source(file.path("calibration", "v02", "pilot_eval.R"))
torch_set_num_threads(1L)
## torch >= 0.14 has the faster ignite Adam (FINN's default); the cluster's torch 0.13 falls back to optim_adam
OPTIMIZER <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
DATA_DIR <- normalizePath(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))[
  dir.exists(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia")))][1])
LABEL <- "conus_C0_v01"; OUT <- file.path("calibration", "v02", "pilots", LABEL); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  prep_version = "conus-annual-v1", train_set = "train", holdout_set = "holdout",
  n_pilot_train = 300L, seed_train = 7L, n_pilot_holdout = 200L, seed_holdout = 8L,   # the fixed pilot subset
  patch_size_ha = 0.01681, n_patches = 4L, recruits_dbh = 12.9,
  env_formula = ~ tavg_c + prec_mm + cwd_z, hidden = c(20L, 20L), dropout = 0.1,
  reg_K_init = 50, reg_floor = 1e-3, adult_x_init = 20, adult_A_init = 15,      # v0.1 settings
  epochs = 250L, batchsize = 100L, lr = 0.01, lr_scheduler = "none",
  weights = "auto", weight_scale = "auto",
  loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  loss_aggregation = list(growth = "tree", mortality = "tree"),
  bare_years = 300L, bare_patches = 10L, seed = 1L
)

## the fixed pilot subset, renumbered 1..n
subset_set <- function(d, n, seed) {
  set.seed(seed); ids <- sort(sample(unique(d$obs$siteID), n)); map <- data.table(siteID = ids, new = seq_along(ids))
  ren <- function(x) { x <- merge(x[siteID %in% ids], map, by = "siteID"); x[, siteID := new][, new := NULL]; setcolorder(x, "siteID"); x[] }
  list(obs = ren(d$obs), env = ren(d$env), init = ren(d$init), trees = ren(d$trees), species = d$species)
}
train <- subset_set(pe_read_set(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)), PARAMS$n_pilot_train, PARAMS$seed_train)
hold  <- subset_set(pe_read_set(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$holdout_set)), PARAMS$n_pilot_holdout, PARAMS$seed_holdout)
n_species <- nrow(train$species)
message(sprintf("pilot subset: %d train plots (%d tree intervals), %d holdout plots", uniqueN(train$obs$siteID), nrow(train$trees), uniqueN(hold$obs$siteID)))

FINN.seed(PARAMS$seed)
model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(PARAMS$env_formula, FINN::growth, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(PARAMS$env_formula, FINN::mortality, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(PARAMS$env_formula, FINN::regeneration_adult, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       custom_parameters = list(reg_logK = rep(log(PARAMS$reg_K_init), n_species),
                                                                reg_adult_logx = rep(log(PARAMS$adult_x_init), n_species),
                                                                reg_adult_logA = rep(log(PARAMS$adult_A_init), n_species))))
t0 <- Sys.time()
fit(model, env = train$env, data = train$obs, tree_data = train$trees, init_cohort = makeInitCohorts(train$init, Nspecies = n_species),
    epochs = PARAMS$epochs, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler,
    weights = PARAMS$weights, weight_scale = PARAMS$weight_scale, loss_family = PARAMS$loss_family, loss_aggregation = PARAMS$loss_aggregation,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, env_autoscale = TRUE, device = "cpu",
    optimizer = OPTIMIZER, folder = NULL, checkpoints = 10L, plot_progress = FALSE)
minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
lh <- as.data.table(do.call(rbind, model$history)); setnames(lh, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(lh))]); lh[, epoch := .I]
fwrite(lh, file.path(OUT, "loss_history.csv"))
torch_save(lapply(model$parameters, function(p) p$detach()$cpu()$clone()), file.path(DATA_DIR, "fits", paste0("pilot_", LABEL, ".pt")))

se <- train$env[, .(tavg_c = mean(tavg_c), prec_mm = mean(prec_mm)), by = siteID][order(tavg_c)]
bare <- se[round(quantile(seq_len(.N), c(0.1, 0.5, 0.9))), .(siteID = .I, tavg_c, prec_mm, cwd_z = 0)]
res <- pilot_evaluate(model, hold, LABEL, OUT, n_patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, step_years = 1L,
                      bare_env = bare, bare_years = PARAMS$bare_years, bare_patches = PARAMS$bare_patches, batchsize = 100L)
K <- exp(as.numeric(as_array(model$reg_logK))); theta <- as.numeric(as_array(model$par_theta_recruits))
res$summary[, `:=`(minutes = round(minutes, 1), epochs = PARAMS$epochs, K_median = median(K), K_max = max(K), theta_median = median(theta),
                   loss_last = paste(sprintf("%s %.3f", names(lh)[1:6], as.numeric(lh[.N, 1:6])), collapse = " "))]
fwrite(res$summary, file.path(OUT, "summary.csv"))
message(sprintf("%s done in %.0f min; K median %.1f max %.0f; theta median %.3f", LABEL, minutes, median(K), max(K), median(theta)))
