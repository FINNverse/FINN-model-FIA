#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/pilot_oregon_O1_cap_cosine.R
##
## WHAT   Pilot O1 = O0 with two settings changed, each tied to an audit finding:
##        1. BOUNDED RECRUITMENT CAP (audit F/G: the v0.1 Beverton-Holt K rose
##           to 1e3-1e5 recruits per ha per decade for 8 species and never
##           bound; bare ground ran to 118-415 m2/ha). FINN's custom parameters
##           can be fixed but not bounded, so K is FIXED per species at
##           2 x the 95th percentile of the observed recruits per ha per decade
##           of that species on the training plots (floor 20), i.e. a cap the
##           data support; the light and climate response below the cap stay
##           free.
##        2. COSINE LEARNING-RATE DECAY (audit F: the v0.1 losses jumped in the
##           last fifth, dbh 0.29 -> 0.97, and no parameter block had settled),
##           lr 0.01 annealed to 0 over the run.
##        Data (oregon-10yr-v2, 18 species), losses, epochs as O0.
## RUN    Rscript calibration/v02/pilot_oregon_O1_cap_cosine.R   (repo root, ~15 min, one thread)
## OUTPUT calibration/v02/pilots/oregon_O1_cap_cosine/
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
source(file.path("calibration", "v02", "pilot_eval.R"))
torch_set_num_threads(1L)
## torch >= 0.14 has the faster ignite Adam (FINN's default); the cluster's torch 0.13 falls back to optim_adam
OPTIMIZER <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
DATA_DIR <- normalizePath(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))[
  dir.exists(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia")))][1])
LABEL <- "oregon_O1_cap_cosine"; OUT <- file.path("calibration", "v02", "pilots", LABEL); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  prep_version = "oregon-10yr-v2", train_set = "train_blocked", holdout_set = "holdout_blocked",   # v0.2 data (as every pilot)
  patch_size_ha = 0.01681, n_patches = 4L, recruits_dbh = 12.9,
  env_formula = ~ temp + prec, hidden = NULL,
  cap_quantile = 0.95, cap_factor = 2, cap_floor = 20,   # K fixed at 2 x q95 of the observed recruits per ha per decade, floor 20
  reg_floor = 0.01,
  epochs = 300L, batchsize = 50L, lr = 0.01, lr_scheduler = "cosine",   # cosine anneal over the run
  weights = "auto", weight_scale = "auto",
  loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  bare_steps = 40L, bare_patches = 10L, seed = 42L
)

train <- pe_read_set(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set))
hold  <- pe_read_set(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$holdout_set))
n_species <- nrow(train$species)
K_cap <- pe_reg_cap(train$obs, train$species, q = PARAMS$cap_quantile, factor = PARAMS$cap_factor, floor_cap = PARAMS$cap_floor, per_step = 1)
message("K cap per species (per ha per decade): ", paste(sprintf("%s %.0f", train$species$species_name, K_cap), collapse = ", "))

FINN.seed(PARAMS$seed)
model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(PARAMS$env_formula, FINN::growth, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(PARAMS$env_formula, FINN::mortality, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(PARAMS$env_formula, FINN::regeneration_saturation, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       custom_parameters = list(reg_logK = list(init = log(K_cap), trainable = FALSE))))
t0 <- Sys.time()
fit(model, env = train$env, data = train$obs, init_cohort = makeInitCohorts(train$init, Nspecies = n_species),
    epochs = PARAMS$epochs, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler,
    lr_scheduler_params = list(T_max = PARAMS$epochs),
    weights = PARAMS$weights, weight_scale = PARAMS$weight_scale, loss_family = PARAMS$loss_family,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, env_autoscale = TRUE, device = "cpu",
    optimizer = OPTIMIZER, folder = NULL, checkpoints = 10L, plot_progress = FALSE)
minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
lh <- as.data.table(do.call(rbind, model$history)); setnames(lh, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(lh))]); lh[, epoch := .I]
fwrite(lh, file.path(OUT, "loss_history.csv"))
torch_save(lapply(model$parameters, function(p) p$detach()$cpu()$clone()), file.path(DATA_DIR, "fits", paste0("pilot_", LABEL, ".pt")))

se <- train$env[, .(temp = mean(temp), prec = mean(prec)), by = siteID][order(temp)]
bare <- se[round(quantile(seq_len(.N), c(0.1, 0.5, 0.9))), .(siteID = .I, temp, prec)]
res <- pilot_evaluate(model, hold, LABEL, OUT, n_patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, step_years = 10L,
                      bare_env = bare, bare_years = PARAMS$bare_steps, bare_patches = PARAMS$bare_patches)
pilot_stratify(res$predictions, hold$sites, OUT)
theta <- as.numeric(as_array(model$par_theta_recruits))
res$summary[, `:=`(minutes = round(minutes, 1), epochs = PARAMS$epochs, K_median = median(K_cap), K_max = max(K_cap), theta_median = median(theta),
                   loss_last = paste(sprintf("%s %.3f", names(lh)[1:6], as.numeric(lh[.N, 1:6])), collapse = " "))]
fwrite(res$summary, file.path(OUT, "summary.csv"))
message(sprintf("%s done in %.0f min; theta median %.3f", LABEL, minutes, median(theta)))
