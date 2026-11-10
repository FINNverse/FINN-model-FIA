#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/pilot_oregon_O3_v2_staged.R
##
## WHAT   Pilot O3 = O2 (data v2 without `other`, fixed data-derived
##        recruitment cap, frozen dispersion, cosine lr) fitted in STAGES
##        instead of one joint run (brief trap 9, loss competition: in the
##        v0.1 joint fit the state losses pulled the regeneration parameters
##        -- `other` and the firs received the recruits -- while the mortality
##        loss never left its intercept baseline):
##          stage A  growth, mortality, competition and their climate
##                   responses on the state losses + growth + mortality
##                   (regeneration weight 0, regeneration parameters frozen:
##                   its climate layer, light response, dispersion, cap);
##          stage B  regeneration alone (weights 0,0,0,0,0,1) with everything
##                   else frozen, so the recruitment process is fitted on its
##                   own likelihood against the observed stand drift of stage A.
##        Both stages anneal their lr with a cosine schedule. 200 + 100 epochs.
## RUN    Rscript calibration/v02/pilot_oregon_O3_v2_staged.R   (repo root, ~15 min, one thread)
## OUTPUT calibration/v02/pilots/oregon_O3_v2_staged/
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
source(file.path("calibration", "v02", "pilot_eval.R"))
source(file.path("R", "definition.R"))
torch_set_num_threads(1L)
## torch >= 0.14 has the faster ignite Adam (FINN's default); the cluster's torch 0.13 falls back to optim_adam
OPTIMIZER <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
DATA_DIR <- normalizePath(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))[
  dir.exists(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia")))][1])
LABEL <- "oregon_O3_v2_staged"; OUT <- file.path("calibration", "v02", "pilots", LABEL); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  prep_version = "oregon-10yr-v2", train_set = "train_blocked", holdout_set = "holdout_blocked",
  patch_size_ha = 0.01681, n_patches = 4L, recruits_dbh = 12.9,
  env_formula = ~ temp + prec, hidden = NULL,
  cap_quantile = 0.95, cap_factor = 2, cap_floor = 20,
  reg_floor = 0.01, theta_source = "oregon", theta_bounds = c(0.05, 2),
  epochs_A = 200L, epochs_B = 100L, batchsize = 50L, lr = 0.01, lr_scheduler = "cosine",
  weights_A = c(1, 1, 1, 1, 1, 0), weights_B = c(0, 0, 0, 0, 0, 1), weight_scale = "auto",
  loss_family = c(dbh = "mse", ba = "mse", trees = "poisson", growth = "mse", mortality = "binomial", regeneration = "nbinom"),
  bare_steps = 40L, bare_patches = 10L, seed = 42L
)

train <- pe_read_set(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$train_set))
hold  <- pe_read_set(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$holdout_set))
n_species <- nrow(train$species)
K_cap <- pe_reg_cap(train$obs, train$species, q = PARAMS$cap_quantile, factor = PARAMS$cap_factor, floor_cap = PARAMS$cap_floor, per_step = 1)
m01 <- build_model(PARAMS$theta_source, torch_load(file.path(DATA_DIR, "fits", PARAMS$theta_source, "parameters.pt")),
                   env_scaling = as.data.frame(fread(file.path("calibration", PARAMS$theta_source, "env_scaling.csv"))))
theta0 <- pmin(pmax(as.numeric(as_array(m01$par_theta_recruits))[seq_len(n_species)], PARAMS$theta_bounds[1]), PARAMS$theta_bounds[2])

FINN.seed(PARAMS$seed)
model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(PARAMS$env_formula, FINN::growth, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(PARAMS$env_formula, FINN::mortality, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  regeneration_process = createProcess(PARAMS$env_formula, FINN::regeneration_saturation, hidden = PARAMS$hidden, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       dispersion_parameter = theta0,
                                       custom_parameters = list(reg_logK = list(init = log(K_cap), trainable = FALSE))))
model$par_theta_recruits_raw$requires_grad_(FALSE)
## the climate layers are created lazily: a one-site dummy simulation instantiates them so they can be frozen per stage
dummy <- data.table(siteID = 1L, year = 1L, temp = mean(train$env$temp), prec = mean(train$env$prec))
model$env_scaling <- NULL
invisible(simulateForest(model, env = dummy, init_cohort = NULL, patches = 1L, patch_size = PARAMS$patch_size_ha, device = "cpu"))
model$env_scaling <- NULL                      # fit() learns the scaling from the training env
set_trainable <- function(model, pattern, on) for (nm in names(model$parameters)) if (grepl(pattern, nm)) model$parameters[[nm]]$requires_grad_(on)
reg_pattern <- "^nn_regeneration|^par_regeneration|^par_theta|^reg_logK"
init_cohort <- makeInitCohorts(train$init, Nspecies = n_species)
t0 <- Sys.time()
## stage A: demography without recruitment learning
set_trainable(model, reg_pattern, FALSE)
fit(model, env = train$env, data = train$obs, init_cohort = init_cohort,
    epochs = PARAMS$epochs_A, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler,
    lr_scheduler_params = list(T_max = PARAMS$epochs_A),
    weights = PARAMS$weights_A, weight_scale = PARAMS$weight_scale, loss_family = PARAMS$loss_family,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, env_autoscale = TRUE, device = "cpu",
    optimizer = OPTIMIZER, folder = NULL, checkpoints = 10L, plot_progress = FALSE)
lhA <- as.data.table(do.call(rbind, model$history)); setnames(lhA, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(lhA))]); lhA[, `:=`(epoch = .I, stage = "A")]
## stage B: recruitment alone, everything else frozen (theta and K stay frozen)
set_trainable(model, "^nn_regeneration|^par_regeneration", TRUE)
set_trainable(model, "^nn_growth|^nn_mortality|^par_growth|^par_mortality|^par_competition", FALSE)
fit(model, env = train$env, data = train$obs, init_cohort = init_cohort,
    epochs = PARAMS$epochs_B, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler,
    lr_scheduler_params = list(T_max = PARAMS$epochs_B),
    weights = PARAMS$weights_B, weight_scale = PARAMS$weight_scale, loss_family = PARAMS$loss_family,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, env_autoscale = TRUE, device = "cpu",
    optimizer = OPTIMIZER, folder = NULL, checkpoints = 10L, plot_progress = FALSE)
minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
lhB <- as.data.table(do.call(rbind, model$history)); setnames(lhB, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(lhB))]); lhB[, `:=`(epoch = .I, stage = "B")]
lh <- rbind(lhA, lhB)
fwrite(lh, file.path(OUT, "loss_history.csv"))
set_trainable(model, ".", TRUE)
torch_save(lapply(model$parameters, function(p) p$detach()$cpu()$clone()), file.path(DATA_DIR, "fits", paste0("pilot_", LABEL, ".pt")))

se <- train$env[, .(temp = mean(temp), prec = mean(prec)), by = siteID][order(temp)]
bare <- se[round(quantile(seq_len(.N), c(0.1, 0.5, 0.9))), .(siteID = .I, temp, prec)]
res <- pilot_evaluate(model, hold, LABEL, OUT, n_patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, step_years = 10L,
                      bare_env = bare, bare_years = PARAMS$bare_steps, bare_patches = PARAMS$bare_patches)
theta <- as.numeric(as_array(model$par_theta_recruits))
res$summary[, `:=`(minutes = round(minutes, 1), epochs = PARAMS$epochs_A + PARAMS$epochs_B, K_median = median(K_cap), K_max = max(K_cap), theta_median = median(theta),
                   loss_last = paste(sprintf("%s %.3f", names(lh)[1:6], as.numeric(lh[.N, 1:6])), collapse = " "))]
fwrite(res$summary, file.path(OUT, "summary.csv"))
message(sprintf("%s done in %.0f min", LABEL, minutes))
