#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/pilot_conus_C3_floor.R
##
## WHAT   Pilot C3 = C2 (data v2, fixed data-derived cap, frozen dispersion,
##        cosine lr, staged demography-then-recruitment, same 300 + 200 pilot
##        plots) with the MORTALITY FLOOR removed: FINN 0.2.0 clamps the
##        simulated annual mortality at >= 1 % (audit finding J); on FIA the
##        observed annual rates of Juniperus osteosperma, Pinus edulis and
##        Quercus alba are 0.3-0.6 %/yr, so the floor binds in the state
##        losses and in every simulation (the per-tree loss is unclamped). The
##        FINN branch `zoo-v02-mortfloor` (FINNverse/FINN#51, PR #52) makes
##        the clamp an argument: `createProcess(mortality_floor = )`. This arm
##        sets mortality_floor = 1e-6 and mortality_draw = "poisson" in the
##        mortality process, nothing else.
## LIBRARY  ~/R/lib-finn-zoo-v02 (the branch build), put FIRST on .libPaths(),
##        before ~/R/lib-finn-0848b60. The script stops if that build lacks the
##        argument. On the cluster the coordinator installs the same build into
##        the library of the same name.
## RUN    Rscript calibration/v02/pilot_conus_C3_floor.R   (repo root, ~40 min, one thread)
## OUTPUT calibration/v02/pilots/conus_C3_floor/
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-zoo-v02", "~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })
if (!"mortality_floor" %in% names(formals(FINN::createProcess)))
  stop("this FINN has no createProcess(mortality_floor=): install the zoo-v02-mortfloor branch build into ~/R/lib-finn-zoo-v02")
message("FINN from ", system.file(package = "FINN"))
source(file.path("calibration", "v02", "pilot_eval.R"))
source(file.path("R", "definition.R"))
torch_set_num_threads(1L)
## torch >= 0.14 has the faster ignite Adam (FINN's default); the cluster's torch 0.13 falls back to optim_adam
OPTIMIZER <- if (exists("optim_ignite_adam", envir = asNamespace("torch"))) torch::optim_ignite_adam else torch::optim_adam
DATA_DIR <- normalizePath(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))[
  dir.exists(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia")))][1])
LABEL <- "conus_C3_floor"; OUT <- file.path("calibration", "v02", "pilots", LABEL); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  prep_version = "conus-annual-v2", train_set = "train", holdout_set = "holdout",
  n_pilot_train = 300L, seed_train = 7L, n_pilot_holdout = 200L, seed_holdout = 8L,
  patch_size_ha = 0.01681, n_patches = 4L, recruits_dbh = 12.9,
  env_formula = ~ tavg_c + prec_mm + cwd_z, hidden = c(20L, 20L), dropout = 0.1,
  cap_quantile = 0.95, cap_factor = 2, cap_floor = 2,
  reg_floor = 1e-3, adult_x_init = 20, adult_A_init = 15,
  theta_source = "conus", theta_bounds = c(0.05, 2),
  mortality_floor = 1e-6, mortality_draw = "poisson",       # THE ONE CHANGE against C2 (FINN#51)
  epochs_A = 170L, epochs_B = 80L, batchsize = 100L, lr = 0.01, lr_scheduler = "cosine",
  weights_A = c(1, 1, 1, 1, 1, 0), weights_B = c(0, 0, 0, 0, 0, 1), weight_scale = "auto",
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

FINN.seed(PARAMS$seed)
model <- finn(
  N_species = n_species, recruits_dbh = PARAMS$recruits_dbh, reg_floor = PARAMS$reg_floor,
  competition_process  = createProcess(~0, FINN::competition, optimizeSpecies = TRUE),
  growth_process       = createProcess(PARAMS$env_formula, FINN::growth, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE),
  mortality_process    = createProcess(PARAMS$env_formula, FINN::mortality, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       mortality_floor = PARAMS$mortality_floor, mortality_draw = PARAMS$mortality_draw),
  regeneration_process = createProcess(PARAMS$env_formula, FINN::regeneration_adult, hidden = PARAMS$hidden, dropout = PARAMS$dropout, optimizeSpecies = TRUE, optimizeEnv = TRUE,
                                       dispersion_parameter = theta0,
                                       custom_parameters = list(reg_logK = list(init = log(K_cap), trainable = FALSE),
                                                                reg_adult_logx = rep(log(PARAMS$adult_x_init), n_species),
                                                                reg_adult_logA = rep(log(PARAMS$adult_A_init), n_species))))
model$par_theta_recruits_raw$requires_grad_(FALSE)
dummy <- data.table(siteID = 1L, year = 1L, tavg_c = mean(train$env$tavg_c), prec_mm = mean(train$env$prec_mm), cwd_z = 0)
invisible(simulateForest(model, env = dummy, init_cohort = NULL, patches = 1L, patch_size = PARAMS$patch_size_ha, device = "cpu"))
model$env_scaling <- NULL
set_trainable <- function(model, pattern, on) for (nm in names(model$parameters)) if (grepl(pattern, nm)) model$parameters[[nm]]$requires_grad_(on)
init_cohort <- makeInitCohorts(train$init, Nspecies = n_species)
t0 <- Sys.time()
set_trainable(model, "^nn_regeneration|^par_regeneration|^par_theta|^reg_", FALSE)
fit(model, env = train$env, data = train$obs, tree_data = train$trees, init_cohort = init_cohort,
    epochs = PARAMS$epochs_A, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler, lr_scheduler_params = list(T_max = PARAMS$epochs_A),
    weights = PARAMS$weights_A, weight_scale = PARAMS$weight_scale, loss_family = PARAMS$loss_family, loss_aggregation = PARAMS$loss_aggregation,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, env_autoscale = TRUE, device = "cpu", optimizer = OPTIMIZER, folder = NULL, checkpoints = 10L, plot_progress = FALSE)
lhA <- as.data.table(do.call(rbind, model$history)); setnames(lhA, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(lhA))]); lhA[, `:=`(epoch = .I, stage = "A")]
set_trainable(model, "^nn_regeneration|^par_regeneration|^reg_adult", TRUE)
set_trainable(model, "^nn_growth|^nn_mortality|^par_growth|^par_mortality|^par_competition", FALSE)
fit(model, env = train$env, data = train$obs, tree_data = train$trees, init_cohort = init_cohort,
    epochs = PARAMS$epochs_B, batchsize = PARAMS$batchsize, lr = PARAMS$lr, lr_scheduler = PARAMS$lr_scheduler, lr_scheduler_params = list(T_max = PARAMS$epochs_B),
    weights = PARAMS$weights_B, weight_scale = PARAMS$weight_scale, loss_family = PARAMS$loss_family, loss_aggregation = PARAMS$loss_aggregation,
    patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, env_autoscale = TRUE, device = "cpu", optimizer = OPTIMIZER, folder = NULL, checkpoints = 10L, plot_progress = FALSE)
minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
lhB <- as.data.table(do.call(rbind, model$history)); setnames(lhB, c("dbh", "ba", "trees", "growth", "mort", "reg", "total")[seq_len(ncol(lhB))]); lhB[, `:=`(epoch = .I, stage = "B")]
lh <- rbind(lhA, lhB); fwrite(lh, file.path(OUT, "loss_history.csv"))
set_trainable(model, ".", TRUE)
p <- lapply(model$parameters, function(x) x$detach()$cpu()$clone()); p$reg_logK <- model$reg_logK$detach()$cpu()$clone()
torch_save(p, file.path(DATA_DIR, "fits", paste0("pilot_", LABEL, ".pt")))

se <- train$env[, .(tavg_c = mean(tavg_c), prec_mm = mean(prec_mm)), by = siteID][order(tavg_c)]
bare <- se[round(quantile(seq_len(.N), c(0.1, 0.5, 0.9))), .(siteID = .I, tavg_c, prec_mm, cwd_z = 0)]
res <- pilot_evaluate(model, hold, LABEL, OUT, n_patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, step_years = 1L,
                      bare_env = bare, bare_years = PARAMS$bare_years, bare_patches = PARAMS$bare_patches, batchsize = 100L)
theta <- as.numeric(as_array(model$par_theta_recruits))
res$summary[, `:=`(minutes = round(minutes, 1), epochs = PARAMS$epochs_A + PARAMS$epochs_B, K_median = median(K_cap), K_max = max(K_cap), theta_median = median(theta),
                   loss_last = paste(sprintf("%s %.3f", names(lh)[1:6], as.numeric(lh[.N, 1:6])), collapse = " "),
                   finn_library = system.file(package = "FINN"), mortality_floor = PARAMS$mortality_floor)]
fwrite(res$summary, file.path(OUT, "summary.csv"))
message(sprintf("%s done in %.0f min", LABEL, minutes))
