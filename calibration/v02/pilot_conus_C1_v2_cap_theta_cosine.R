#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/pilot_conus_C1_v2_cap_theta_cosine.R
##
## WHAT   Pilot C1 = C0 (same 300 + 200 pilot plots, same model) with the
##        settings changed that the audit motivates:
##        1. DATA v2 (conus-annual-v2, calibration/v02/prepare_conus_v2*.R): the
##           recruitment target is the complete ingrowth into 12.7 cm (audit B:
##           v0.1 missed the 8 % tracked from the microplot and double-counted
##           re-recorded trees; train +3.6 %, holdout -4.4 %).
##        2. FIXED RECRUITMENT CAP (brief trap 8): K per species fixed at 2 x the
##           95th percentile of the observed annual recruits per ha of that
##           species (floor 2 per ha per year) instead of a trainable K
##           (audit F: K was still moving at the end of v0.1, 0.85 of its
##           early change).
##        3. FROZEN DISPERSION (brief trap 7): theta per species from the v0.1
##           production fit, clipped to [0.05, 2], not trained.
##        4. COSINE lr decay (audit F: the regeneration network's change per
##           snapshot was still GROWING at the end of v0.1, ratio 2.5-7).
## RUN    Rscript calibration/v02/pilot_conus_C1_v2_cap_theta_cosine.R   (repo root, ~40 min, one thread)
## OUTPUT calibration/v02/pilots/conus_C1_v2_cap_theta_cosine/
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
source(file.path("calibration", "v02", "pilot_eval.R"))
source(file.path("R", "definition.R"))
torch_set_num_threads(1L)
DATA_DIR <- normalizePath(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))[
  dir.exists(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia")))][1])
LABEL <- "conus_C1_v2_cap_theta_cosine"; OUT <- file.path("calibration", "v02", "pilots", LABEL); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  prep_version = "conus-annual-v2", train_set = "train", holdout_set = "holdout",
  n_pilot_train = 300L, seed_train = 7L, n_pilot_holdout = 200L, seed_holdout = 8L,   # the same pilot plots as C0 (same seeds, same site ids)
  patch_size_ha = 0.01681, n_patches = 4L, recruits_dbh = 12.9,
  env_formula = ~ tavg_c + prec_mm + cwd_z, hidden = c(20L, 20L), dropout = 0.1,
  cap_quantile = 0.95, cap_factor = 2, cap_floor = 2,      # per ha per YEAR: observed reg / period_length
  reg_floor = 1e-3, adult_x_init = 20, adult_A_init = 15,
  theta_source = "conus", theta_bounds = c(0.05, 2),
  epochs = 250L, batchsize = 100L, lr = 0.01, lr_scheduler = "cosine",
  weights = "auto", weight_scale = "auto",
  loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  loss_aggregation = list(growth = "tree", mortality = "tree"),
  bare_years = 300L, bare_patches = 10L, seed = 1L
)
subset_set <- function(d, n, seed) {
  set.seed(seed); ids <- sort(sample(unique(d$obs$siteID), n)); map <- data.table(siteID = ids, new = seq_along(ids))
  ren <- function(x) { x <- merge(x[siteID %in% ids], map, by = "siteID"); x[, siteID := new][, new := NULL]; setcolorder(x, "siteID"); x[] }
  list(obs = ren(d$obs), env = ren(d$env), init = ren(d$init), trees = ren(d$trees), species = d$species)
}
train <- subset_set(pe_read_set(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set)), PARAMS$n_pilot_train, PARAMS$seed_train)
hold  <- subset_set(pe_read_set(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$holdout_set)), PARAMS$n_pilot_holdout, PARAMS$seed_holdout)
n_species <- nrow(train$species)
obs_annual <- copy(train$obs)[, reg := reg / period_length]
K_cap <- pe_reg_cap(obs_annual, train$species, q = PARAMS$cap_quantile, factor = PARAMS$cap_factor, floor_cap = PARAMS$cap_floor, per_step = 1)
m01 <- build_model(PARAMS$theta_source, torch_load(file.path(DATA_DIR, "fits", PARAMS$theta_source, "parameters.pt")),
                   env_scaling = as.data.frame(fread(file.path("calibration", PARAMS$theta_source, "env_scaling.csv"))))
theta0 <- pmin(pmax(as.numeric(as_array(m01$par_theta_recruits)), PARAMS$theta_bounds[1]), PARAMS$theta_bounds[2])
message("K cap (per ha per year): ", paste(sprintf("%s %.1f", train$species$species_name, K_cap), collapse = ", "))

FINN.seed(PARAMS$seed)
model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(PARAMS$env_formula, FINN::growth, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(PARAMS$env_formula, FINN::mortality, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(PARAMS$env_formula, FINN::regeneration_adult, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       dispersion_parameter = theta0,
                                       custom_parameters = list(reg_logK = list(init = log(K_cap), trainable = FALSE),
                                                                reg_adult_logx = rep(log(PARAMS$adult_x_init), n_species),
                                                                reg_adult_logA = rep(log(PARAMS$adult_A_init), n_species))))
model$par_theta_recruits_raw$requires_grad_(FALSE)
t0 <- Sys.time()
fit(model, env = train$env, data = train$obs, tree_data = train$trees, init_cohort = makeInitCohorts(train$init, Nspecies = n_species),
    epochs = PARAMS$epochs, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler,
    lr_scheduler_params = list(T_max = PARAMS$epochs),
    weights = PARAMS$weights, weight_scale = PARAMS$weight_scale, loss_family = PARAMS$loss_family, loss_aggregation = PARAMS$loss_aggregation,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, env_autoscale = TRUE, device = "cpu",
    folder = NULL, checkpoints = 10L, plot_progress = FALSE)
minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
lh <- as.data.table(do.call(rbind, model$history)); setnames(lh, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(lh))]); lh[, epoch := .I]
fwrite(lh, file.path(OUT, "loss_history.csv"))
torch_save(lapply(model$parameters, function(p) p$detach()$cpu()$clone()), file.path(DATA_DIR, "fits", paste0("pilot_", LABEL, ".pt")))

se <- train$env[, .(tavg_c = mean(tavg_c), prec_mm = mean(prec_mm)), by = siteID][order(tavg_c)]
bare <- se[round(quantile(seq_len(.N), c(0.1, 0.5, 0.9))), .(siteID = .I, tavg_c, prec_mm, cwd_z = 0)]
res <- pilot_evaluate(model, hold, LABEL, OUT, n_patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, step_years = 1L,
                      bare_env = bare, bare_years = PARAMS$bare_years, bare_patches = PARAMS$bare_patches, batchsize = 100L)
theta <- as.numeric(as_array(model$par_theta_recruits))
res$summary[, `:=`(minutes = round(minutes, 1), epochs = PARAMS$epochs, K_median = median(K_cap), K_max = max(K_cap), theta_median = median(theta),
                   loss_last = paste(sprintf("%s %.3f", names(lh)[1:6], as.numeric(lh[.N, 1:6])), collapse = " "))]
fwrite(res$summary, file.path(OUT, "summary.csv"))
message(sprintf("%s done in %.0f min", LABEL, minutes))
