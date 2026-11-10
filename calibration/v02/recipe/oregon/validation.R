#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-zoo | FINN-model-FIA/calibration/v02/recipe/oregon/validation.R
##
## WHAT   The validation report of the zoo v0.2 recipe model of the Oregon
##        variant (fit_stageA.R + fit_stageB.R, fits/oregon-10yr-v2-recipe/model.pt).
##        Runs on the saved model, so it can be redone after a plot change without
##        a refit. Two sets, named explicitly:
##          TEST    = the plots the model was fitted on (prepared/oregon-10yr-v2/train_blocked, all 450),
##          HOLDOUT = plots the model never saw, in 1-degree blocks that hold no
##                    training plot (prepared/oregon-10yr-v2/holdout_blocked, 111).
##        The model has no coordinates as predictors, so the holdout is predicted
##        exactly like the test set.
##        For every plot the 2001-03 stand is simulated forward over its two
##        10-year steps and each step is scored the way fit() scores it: states at
##        the step's census, growth and stem loss as the rate per step, recruitment
##        per ha per step. Metrics per response and set (species rows and stand
##        totals) with the PERSISTENCE null (end = the state at the start of the
##        step: the first inventory, or the previous census) and the change over
##        the step beside them, rate calibration ratios (predicted / observed
##        recruits, deaths, growth), holdout errors by the END-census out-of-S
##        share, mortality AUC per site weighted by trees at risk (no tree table in
##        this variant, so no per-tree AUC), a per-species Spearman table, and the
##        figures the zoo card shows: loss curves of both stages, observed vs
##        predicted, site-level ROC, map, bare-ground succession across the climate
##        space (50 steps = 500 years), parameter trajectories of both stages.
##        Copied from FINN-model-IFN-Spain/calibration/v02/recipe/validation.R
##        (cb41778) and adapted to the decade step, two steps per plot and the
##        absence of tree observations.
## INPUT  <FIT_DIR>/{model.pt, model_stageA.pt, parameters.pt, loss_history_stageA.csv,
##        loss_history_stageB.csv, fit_summary_stage{A,B}.md}; prepared/oregon-10yr-v2/{train_blocked,holdout_blocked}
## OUTPUT <FIT_DIR>/report.md, <FIT_DIR>/validation/*.csv, *.png;
##        copies of metrics.csv, species.csv and the card figures into the zoo
##        folder when it exists (ZOO_DIR below)
## RUN    Rscript calibration/v02/recipe/oregon/validation.R        (from the repo root, 4 threads)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
if (dir.exists("~/R/lib-finn-0848b60")) .libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({
  library(FINN); library(data.table); library(torch); library(ggplot2)
})
stopifnot(packageVersion("FINN") == "0.2.0")

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "fia")
FIT_DIR  <- file.path(DATA_DIR, "fits", "oregon-10yr-v2-recipe")
VAL      <- file.path(FIT_DIR, "validation")
ZOO_DIR  <- file.path("..", "..", "FINN-zoo-worktrees", "fia", "FIA-oregon")   # model-card folder (optional)
dir.create(VAL, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version  = "oregon-10yr-v2",
  sets          = list(test = "train_blocked", holdout = "holdout_blocked"),
  n_test_sites  = 450L,                    # subset of the training plots simulated for the report (all of them)
  subset_seed   = 7L,
  patch_size_ha = 0.01681,                 # as fitted
  n_patches     = 4L,                      # as fitted: the four subplots
  env_vars      = c("temp", "prec"),       # as fitted (static WorldClim normals)
  step_years    = 10L,                     # one simulation step = one remeasurement interval
  sim_seed      = 1L,
  ## bare-ground succession: plots nearest to the 6 corners / centre of the
  ## (temp, prec) space of the training set, simulated for eq_steps steps
  eq_steps      = 50L,                     # 50 steps = 500 years
  eq_patches    = 20L,
  eq_quantiles  = list(c(0.1, 0.1), c(0.1, 0.9), c(0.5, 0.5), c(0.9, 0.1), c(0.9, 0.9), c(0.3, 0.7)),
  n_top_species = 5L,                      # species coloured in the succession figure
  threads       = 4L
)
torch_set_num_threads(PARAMS$threads)
PS <- PARAMS$patch_size_ha
ba_tree <- function(dbh) pi * (dbh / 200)^2

## ---- 0. model and data -------------------------------------------------------
model <- torch_load(file.path(FIT_DIR, "model.pt"))
model$eval()
stageA <- torch_load(file.path(FIT_DIR, "model_stageA.pt"))
## the parameter dump must name exactly the tensors of the saved model (the package's
## build_model() rebuilds the model from it; its structure is the recipe's)
parameters <- torch_load(file.path(FIT_DIR, "parameters.pt"))
sd_model <- model$state_dict()
stopifnot(setequal(names(parameters$state_dict), names(sd_model)),
          all(vapply(names(sd_model), function(n) identical(as.integer(parameters$state_dict[[n]]$shape), as.integer(sd_model[[n]]$shape)), logical(1))),
          max(vapply(names(sd_model), function(n) max(abs(as.numeric(parameters$state_dict[[n]]) - as.numeric(sd_model[[n]]))), numeric(1))) == 0)
message("parameters.pt: every tensor of model.pt present, same shapes and values")
loss_hist <- rbind(fread(file.path(FIT_DIR, "loss_history_stageA.csv"))[, stage := "A: demography"],
                   fread(file.path(FIT_DIR, "loss_history_stageB.csv"))[, stage := "B: regeneration"], fill = TRUE)
species_dt <- fread(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$sets$test, "species_dt.csv"))
n_species <- nrow(species_dt)

read_set <- function(set, ids = NULL) {
  d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, set)
  out <- list(obs = fread(file.path(d, "obs_dt.csv")), env = fread(file.path(d, "env_dt.csv"))[, c("siteID", "year", PARAMS$env_vars), with = FALSE],
              init = fread(file.path(d, "init_cohort_dt.csv")), sites = fread(file.path(d, "siteID_dt.csv")))
  if (!"period_length" %in% names(out$obs)) out$obs[, period_length := 1L]   # one step per observation
  if (!is.null(ids)) {
    ## FINN wants siteIDs 1..n: renumber the subset consistently in every table
    key <- data.table(siteID = ids, new = seq_along(ids))
    out <- lapply(out, function(x) {
      x <- merge(x[siteID %in% ids], key, by = "siteID")
      x[, siteID := new][, new := NULL]
      setorder(x, siteID)
    })
  }
  ## the END-census out-of-S share bins of the agreed rare-species rule
  out$sites[, oos_bin_end := cut(share_out_end, c(-1, 0.01, 0.05, 2), labels = c("< 1 %", "1-5 %", "> 5 %"))]
  out
}
set.seed(PARAMS$subset_seed)
n_train_all <- nrow(fread(file.path(DATA_DIR, "prepared", PARAMS$prep_version, PARAMS$sets$test, "siteID_dt.csv")))
test_ids <- sort(sample(seq_len(n_train_all), min(PARAMS$n_test_sites, n_train_all)))
sets <- list(test = read_set(PARAMS$sets$test, test_ids), holdout = read_set(PARAMS$sets$holdout))

## ---- 1. simulate every plot over its steps and aggregate like fit() ----------
## (simulateForest() labels the step driven by env row k as year k; the observation at
## y_end is driven by rows y_end - p + 1 .. y_end, here the one row y_end)
intervals <- function(d) unique(d$obs[, .(siteID, y_end = year, y0 = year - period_length, period_length)])[order(siteID, y_end)]
simulate_set <- function(d) {
  ic <- makeInitCohorts(d$init, Nspecies = n_species)
  FINN.seed(PARAMS$sim_seed)
  sim <- simulateForest(model, env = d$env, patches = PARAMS$n_patches, patch_size = PS,
                        init_cohort = ic, device = "cpu")$long$site
  for (cc in c("siteID", "year", "species")) sim[, (cc) := as.integer(as.character(get(cc)))]
  w <- dcast(sim, siteID + year + species ~ variable, value.var = "value")
  iv <- intervals(d)
  W <- iv[, .(year = seq.int(y0 + 1L, y_end)), by = .(siteID, y_end)]
  w <- merge(w, W, by = c("siteID", "year"), allow.cartesian = TRUE)
  w[, at_end := year == y_end]   # marked before grouping: inside `by = .(year = y_end)` the name `year` is the group's scalar
  out <- w[, .(ba = ba[at_end], trees = trees[at_end], dbh = dbh[at_end],
               growth = mean(growth), mort = mean(mort), reg = sum(r_mean_ha)), by = .(siteID, year = y_end, species)]
  stopifnot(nrow(out) == nrow(iv) * n_species, out[, .N, by = .(siteID, year, species)][, all(N == 1)])
  out
}
pred <- lapply(sets, simulate_set)
## the state at the START of every step (persistence null): the first inventory for
## step 1, the previous census for step 2 (per patch, like obs)
start_state <- function(d) {
  iv <- intervals(d)
  init <- if ("living" %in% names(d$init)) d$init[living == TRUE] else d$init
  s0 <- init[trees > 0, .(ba = sum(trees * ba_tree(dbh)) / PARAMS$n_patches, trees = sum(trees) / PARAMS$n_patches, dbh = weighted.mean(dbh, trees)), by = .(siteID, species)]
  s0 <- merge(s0, iv[, .(year = min(y_end)), by = siteID], by = "siteID")
  s1 <- merge(iv[y0 > 0, .(siteID, y0, year = y_end)], d$obs[, .(siteID, y0 = year, species, ba, trees, dbh)], by = c("siteID", "y0"), allow.cartesian = TRUE)
  rbind(s0[, .(siteID, year, species, ba, trees, dbh)], s1[, .(siteID, year, species, ba, trees, dbh)])
}
start <- lapply(sets, start_state)

## observed vs predicted, long
vars <- c("ba", "trees", "dbh", "growth", "mort", "reg")
OP <- rbindlist(lapply(names(sets), function(s) {
  o <- melt(sets[[s]]$obs[, c("siteID", "year", "species", vars), with = FALSE], id.vars = c("siteID", "year", "species"),
            variable.name = "variable", value.name = "obs")
  p <- melt(pred[[s]][, c("siteID", "year", "species", vars), with = FALSE], id.vars = c("siteID", "year", "species"),
            variable.name = "variable", value.name = "pred")
  merge(o, p, by = c("siteID", "year", "species", "variable"))[, set := s]
}))
OP <- merge(OP, species_dt, by = "species")
## stand totals per step: ba and trees summed, dbh as the stem-weighted mean; persistence = the start state
TOT <- rbindlist(lapply(names(sets), function(s) {
  tot <- function(x) x[, .(ba = sum(ba), trees = sum(trees), dbh = if (sum(trees) > 0) weighted.mean(dbh, trees, na.rm = TRUE) else NA_real_), by = .(siteID, year)]
  o <- tot(sets[[s]]$obs)[, reg := sets[[s]]$obs[, sum(reg), by = .(siteID, year)]$V1]
  p <- tot(pred[[s]])[, reg := pred[[s]][, sum(reg), by = .(siteID, year)]$V1]
  q <- tot(start[[s]])[, reg := NA_real_]
  m <- merge(melt(o, id.vars = c("siteID", "year"), variable.name = "variable", value.name = "obs"),
             melt(p, id.vars = c("siteID", "year"), variable.name = "variable", value.name = "pred"), by = c("siteID", "year", "variable"))
  merge(m, melt(q, id.vars = c("siteID", "year"), variable.name = "variable", value.name = "persist"), by = c("siteID", "year", "variable"), all.x = TRUE)[, set := s]
}))

metrics <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred); o <- obs[ok]; p <- pred[ok]
  if (length(o) < 3) return(list(n = length(o), rmse = NA_real_, r2 = NA_real_, r = NA_real_, bias = NA_real_, obs_mean = NA_real_))
  list(n = length(o), rmse = sqrt(mean((o - p)^2)), r2 = 1 - sum((o - p)^2) / sum((o - mean(o))^2),
       r = suppressWarnings(cor(o, p)), bias = mean(p - o), obs_mean = mean(o))
}
## states per patch -> per ha for the report
OP[variable %in% c("ba", "trees"), `:=`(obs = obs / PS, pred = pred / PS)]
TOT[variable %in% c("ba", "trees"), `:=`(obs = obs / PS, pred = pred / PS, persist = persist / PS)]
M_species <- OP[, metrics(obs, pred), by = .(set, variable)][, level := "species rows (site x step x species)"]
M_total   <- TOT[, metrics(obs, pred), by = .(set, variable)][, level := "stand totals (site x step)"]
M_persist <- TOT[variable %in% c("ba", "trees", "dbh"), metrics(obs, persist), by = .(set, variable)][, level := "stand totals: persistence (end = start of the step)"]
M_change  <- TOT[variable %in% c("ba", "trees"), metrics(obs - persist, pred - persist), by = .(set, variable)][, `:=`(level = "stand totals: change over the step", variable = paste0("d_", variable))]
M <- rbind(M_species, M_total, M_persist, M_change)
spearman <- function(o, p) { ok <- is.finite(o) & is.finite(p); if (sum(ok) < 3) NA_real_ else suppressWarnings(cor(o[ok], p[ok], method = "spearman")) }
M_bysp <- OP[, .(n = sum(is.finite(obs) & is.finite(pred)), spearman = spearman(obs, pred),
                 bias = mean(pred - obs, na.rm = TRUE), obs_mean = mean(obs, na.rm = TRUE)), by = .(set, variable, species_name)]
fwrite(M, file.path(VAL, "metrics.csv")); fwrite(M_bysp, file.path(VAL, "metrics_by_species.csv"))
## holdout errors by the END-census out-of-S share (the agreed rare-species rule's stratified evaluation)
BIN <- merge(TOT[set == "holdout" & variable %in% c("ba", "trees")], sets$holdout$sites[, .(siteID, oos_bin_end)], by = "siteID")
M_bin <- BIN[, c(metrics(obs, pred), list(rmse_persist = sqrt(mean((obs - persist)^2)))), by = .(variable, oos_bin_end)][order(variable, oos_bin_end)]
M_bin[, rmse_ratio_model_persist := rmse / rmse_persist]
fwrite(M_bin, file.path(VAL, "metrics_by_oos_bin.csv"))

## ---- 2. rates: AUC per site, calibration ratios --------------------------------------
auc_weighted <- function(p, pos, neg) {           # ties count half; pos/neg = counts sharing one probability
  o <- order(p); p <- p[o]; pos <- pos[o]; neg <- neg[o]
  cum_neg <- cumsum(neg) - neg                     # negatives with a smaller p
  ties_neg <- neg                                  # negatives at the same p
  (sum(pos * (cum_neg + 0.5 * ties_neg))) / (sum(pos) * sum(neg))
}
## the simulated mortality is the death probability per step (the step IS the interval)
SITE <- rbindlist(lapply(names(sets), function(s) merge(pred[[s]][, .(siteID, year, species, mort_pred = pmin(pmax(mort, 0), 1), growth_pred = growth)],
  sets[[s]]$obs[, .(siteID, year, species, n_at_risk, n_died, growth, growth_n)], by = c("siteID", "year", "species"))[, set := s]))
AUC <- SITE[, {
  m <- .SD[n_at_risk > 0 & is.finite(mort_pred)]
  g <- .SD[growth_n > 0 & is.finite(growth) & is.finite(growth_pred)]
  gm <- metrics(g$growth, g$growth_pred)
  list(site_auc_weighted = m[, auc_weighted(mort_pred, n_died, n_at_risk - n_died)], n_rows = nrow(m), n_at_risk = sum(m$n_at_risk), n_died = sum(m$n_died),
       site_growth_r2 = gm$r2, site_growth_rmse = gm$rmse,
       deaths_ratio = m[, sum(mort_pred * n_at_risk) / sum(n_died)],
       growth_ratio = g[, weighted.mean(growth_pred, growth_n) / weighted.mean(growth, growth_n)],
       reg_ratio = OP[set == .BY$set & variable == "reg", sum(pred, na.rm = TRUE) / sum(obs, na.rm = TRUE)])
}, by = set]
fwrite(SITE, file.path(VAL, "site_rate_predictions.csv"))
fwrite(AUC, file.path(VAL, "mortality_auc.csv"))
## per-species recruitment and death ratios on the holdout
SP <- merge(OP[set == "holdout" & variable == "reg", .(reg_obs = sum(obs, na.rm = TRUE), reg_ratio = round(sum(pred, na.rm = TRUE) / sum(obs, na.rm = TRUE), 2)), by = species_name],
            SITE[set == "holdout" & n_at_risk > 0, .(n_died = sum(n_died), deaths_ratio = round(sum(mort_pred * n_at_risk) / sum(n_died), 2)), by = species][
              , species_name := species_dt$species_name[species]][, .(species_name, n_died, deaths_ratio)], by = "species_name", all.x = TRUE)[order(-reg_obs)]
fwrite(SP, file.path(VAL, "rates_by_species_holdout.csv"))

## ---- 3. figures ------------------------------------------------------------------
theme_set(theme_minimal(base_size = 10))
## 3a loss curves of both stages
LH <- melt(loss_hist[, .(stage, epoch, dbh, ba, trees, growth, mort, reg)], id.vars = c("stage", "epoch"))[is.finite(value) & value != 0]
g1 <- ggplot(LH, aes(epoch, value)) + geom_line(colour = "#1b6ca8") + facet_grid(stage ~ variable, scales = "free") +
  labs(y = "loss (baseline-scaled)", title = "Training loss per response: stage A (demography, observed recruitment) and stage B (regeneration)")
ggsave(file.path(VAL, "fig1_loss.png"), g1, width = 12, height = 4.5, dpi = 130, bg = "white")
## 3b observed vs predicted (species rows + totals)
OPf <- rbind(OP[, .(set, variable, obs, pred, level = "site x step x species")], TOT[, .(set, variable, obs, pred, level = "stand total")])
OPf[, variable := factor(variable, levels = vars, labels = c("basal area (m2/ha)", "stems (/ha)", "mean dbh (cm)", "growth (rel./step)", "stem loss (/step)", "recruits (/ha/step)"))]
g2 <- ggplot(OPf[is.finite(obs) & is.finite(pred)], aes(obs, pred, colour = set)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50") + geom_point(alpha = 0.15, size = 0.6) +
  facet_grid(level ~ variable, scales = "free") + scale_colour_manual(values = c(test = "#1b6ca8", holdout = "#d1495b")) +
  labs(x = "observed", y = "predicted", title = "Observed vs predicted at the census closing each 10-year step", colour = NULL) + theme(legend.position = "bottom")
ggsave(file.path(VAL, "fig2_obs_pred.png"), g2, width = 13, height = 5.5, dpi = 130, bg = "white")
## 3c ROC per site x step x species row, weighted by the trees at risk (no tree table)
roc_dt <- SITE[n_at_risk > 0 & is.finite(mort_pred)][order(set, -mort_pred)][
  , .(tpr = cumsum(n_died) / sum(n_died), fpr = cumsum(n_at_risk - n_died) / sum(n_at_risk - n_died)), by = set]
g3 <- ggplot(roc_dt, aes(fpr, tpr, colour = set)) + geom_abline(slope = 1, intercept = 0, colour = "grey60") + geom_path() +
  scale_colour_manual(values = c(test = "#1b6ca8", holdout = "#d1495b")) + coord_equal() +
  labs(x = "false positive rate", y = "true positive rate", title = "Stem loss per site x step x species: ROC (weighted by trees at risk)",
       subtitle = paste(sprintf("%s AUC %.3f", AUC$set, AUC$site_auc_weighted), collapse = "; "), colour = NULL)
ggsave(file.path(VAL, "fig3_roc.png"), g3, width = 5, height = 5, dpi = 130, bg = "white")
## 3d map: observed and predicted stand basal area at the last census
MAP <- rbindlist(lapply(names(sets), function(s) merge(TOT[set == s & variable == "ba"][, .SD[year == max(year)], by = siteID], sets[[s]]$sites[, .(siteID, lon, lat)], by = "siteID")))
MAPl <- melt(MAP[, .(set, lon, lat, observed = obs, predicted = pred)], id.vars = c("set", "lon", "lat"))
g4 <- ggplot(MAPl, aes(lon, lat, colour = pmin(value, 80))) + geom_point(size = 0.8) + facet_grid(set ~ variable) + coord_quickmap() +
  scale_colour_viridis_c(name = "basal area\n(m2/ha, capped 80)") + labs(title = "Stand basal area at the last census (2021-24)", x = NULL, y = NULL)
ggsave(file.path(VAL, "fig4_map.png"), g4, width = 9, height = 6, dpi = 130, bg = "white")
## 3e bare-ground succession at plots spanning the climate space
train_sites <- merge(sets$test$sites[, .(siteID)], unique(sets$test$env[, c("siteID", PARAMS$env_vars), with = FALSE]), by = "siteID")
pick <- unique(rbindlist(lapply(PARAMS$eq_quantiles, function(q) {
  tq <- quantile(train_sites$temp, q[1]); pq <- quantile(train_sites$prec, q[2])
  z <- train_sites[, (scale(temp)[, 1] - (tq - mean(temp)) / sd(temp))^2 + (scale(prec)[, 1] - (pq - mean(prec)) / sd(prec))^2]
  train_sites[which.min(z)]
})), by = "siteID")
eq_env <- pick[, .(year = seq_len(PARAMS$eq_steps)), by = c("siteID", PARAMS$env_vars)][, c("siteID", "year", PARAMS$env_vars), with = FALSE]
eq_env[, siteID := match(siteID, pick$siteID)]
pick[, siteID := seq_len(.N)]
FINN.seed(PARAMS$sim_seed)
eq <- simulateForest(model, env = eq_env, patches = PARAMS$eq_patches, patch_size = PS, init_cohort = NULL, device = "cpu")$long$site
for (cc in c("siteID", "year", "species")) eq[, (cc) := as.integer(as.character(get(cc)))]
eq <- merge(eq[variable %in% c("ba", "trees")], pick, by = "siteID")
eq <- merge(eq, species_dt, by = "species")
eq[, value := value / PS]
eq[, years := year * PARAMS$step_years]
eq[, site_lab := sprintf("%.1f C, %.0f mm", temp, prec)]
top <- eq[variable == "ba", .(mx = max(value)), by = .(siteID, species_name)][order(siteID, -mx)][, head(.SD, PARAMS$n_top_species), by = siteID]
eq[, shown := paste(siteID, species_name) %in% paste(top$siteID, top$species_name)]
g5 <- ggplot() +
  geom_line(data = eq[variable == "ba" & !shown], aes(years, value, group = species_name), colour = "grey80", linewidth = 0.3) +
  geom_line(data = eq[variable == "ba" & shown], aes(years, value, colour = species_name), linewidth = 0.6) +
  geom_line(data = eq[variable == "ba", .(value = sum(value)), by = .(site_lab, years)], aes(years, value), colour = "black", linewidth = 0.5, linetype = 2) +
  facet_wrap(~ site_lab, ncol = 3) + labs(x = "year", y = "basal area (m2/ha)", colour = NULL,
  title = "Succession from bare ground (dashed: stand total)", subtitle = sprintf("%d patches of %.4f ha per site, %d steps of %d years; the %d species with the highest basal area coloured", PARAMS$eq_patches, PS, PARAMS$eq_steps, PARAMS$step_years, PARAMS$n_top_species)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))
ggsave(file.path(VAL, "fig5_bare_ground.png"), g5, width = 11, height = 7, dpi = 130, bg = "white")
EQ <- eq[variable == "ba", .(ba_total = sum(value)), by = .(siteID, site_lab, year)][
  , .(ba_last10 = mean(ba_total[year > PARAMS$eq_steps - 10]), ba_prev10 = mean(ba_total[year > PARAMS$eq_steps - 20 & year <= PARAMS$eq_steps - 10])), by = .(siteID, site_lab)]
EQ[, rel_change := (ba_last10 - ba_prev10) / ba_prev10]
EQ <- merge(EQ, eq[variable == "trees" & year > PARAMS$eq_steps - 10, .(stems_last10 = sum(value) / 10), by = siteID], by = "siteID")
fwrite(EQ, file.path(VAL, "equilibrium.csv"))
## 3f parameter trajectories of both stages (mechanistic blocks; the networks as their mean change)
traj <- function(ph, n_epochs, stage) {
  if (length(ph) < 2) return(NULL)
  step <- n_epochs / length(ph)
  rbindlist(lapply(seq_along(ph), function(i) {
    s <- ph[[i]]
    rbindlist(lapply(names(s), function(nm) {
      mat <- s[[nm]]
      data.table(stage = stage, epoch = i * step, block = nm, row = as.integer(row(mat)), col = as.integer(col(mat)), value = as.numeric(mat))
    }))
  }))
}
PT <- rbind(traj(stageA$param_history, loss_hist[stage == "A: demography", .N], "A"), traj(model$param_history, loss_hist[stage == "B: regeneration", .N], "B"))
if (!is.null(PT) && nrow(PT)) {
  PT[, kind := fifelse(grepl("^nn_", block), "environment network", "mechanistic")]
  fwrite(PT, file.path(VAL, "parameter_trajectories.csv"))
  setorder(PT, stage, block, row, col, epoch)
  PT[, delta := abs(value - shift(value)), by = .(stage, block, row, col)]
  ## only the blocks a stage trains are shown for it
  moved <- PT[!is.na(delta), .(moved = max(delta) > 0), by = .(stage, block)][moved == TRUE]
  mech <- PT[kind == "mechanistic"][moved, on = .(stage, block)]
  mech[, panel := sprintf("%s: %s [, %d]", stage, sub("_unconstrained", "", block), col)]
  mech[block == "reg_logK", value := exp(value)]
  g6a <- ggplot(mech, aes(epoch, value, group = row)) + geom_line(alpha = 0.35, linewidth = 0.3, colour = "#1b6ca8") +
    facet_wrap(~ panel, scales = "free", ncol = 4) + theme_minimal(base_size = 8) +
    labs(x = "epoch", y = "value", title = "Mechanistic per-species parameters over training (one line per species; A = demography, B = regeneration)")
  conv <- PT[!is.na(delta)][moved, on = .(stage, block)][, .(mean_abs_change = mean(delta)), by = .(stage, epoch, block, kind)]
  g6b <- ggplot(conv, aes(epoch, mean_abs_change, colour = block)) + geom_line(linewidth = 0.4) + scale_y_log10() +
    facet_grid(kind ~ stage, scales = "free") + theme_minimal(base_size = 8) +
    theme(legend.text = element_text(size = 6), legend.key.height = unit(7, "pt")) +
    labs(x = "epoch", y = "mean |change| per snapshot (log)", title = "Are the parameters still moving?")
  ggsave(file.path(VAL, "fig6_parameters.png"), gridExtra::grid.arrange(g6a, g6b, ncol = 1, heights = c(1.3, 1)),
         width = 12, height = 11, dpi = 130, bg = "white")
}

## ---- 4. report -------------------------------------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "-", formatC(x, digits = d, format = "f"))
Mtab <- M[, .(set, level, variable, n, RMSE = fmt(rmse), R2 = fmt(r2), r = fmt(r), bias = fmt(bias), obs_mean = fmt(obs_mean))]
summ <- c(fread(file.path(FIT_DIR, "fit_summary_stageA.md"), sep = NULL, header = FALSE)$V1, "",
          fread(file.path(FIT_DIR, "fit_summary_stageB.md"), sep = NULL, header = FALSE)$V1)
provenance <- data.table(
  script = c("calibration/oregon/01_prepare.R", "calibration/v02/species_tradeoff.R", "calibration/v02/prepare_oregon_v2.R",
             "calibration/v02/recipe/oregon/fit_stageA.R", "calibration/v02/recipe/oregon/fit_stageB.R", "calibration/v02/recipe/oregon/validation.R"),
  output = c("prepared/oregon-10yr-v1 (the Oregon export, FIADB coordinates, WorldClim normals; calibration/PROVENANCE.md)",
             "calibration/v02/species_set_oregon.csv (S = 18, the agreed rare-species rule), species_tradeoff.csv",
             sprintf("prepared/%s/{all,train_blocked,holdout_blocked}, PREPARATION.md, calibration/v02/expected_preparation_oregon_v2.csv", PARAMS$prep_version),
             "fits/oregon-10yr-v2-recipe/{model_stageA.pt, loss_history_stageA.csv, species_parameters_stageA.csv, fit_summary_stageA.md}",
             "fits/oregon-10yr-v2-recipe/{model.pt, parameters.pt, loss_history_stageB.csv, regeneration_parameters.csv, fit_summary_stageB.md}",
             "fits/oregon-10yr-v2-recipe/{report.md, validation/}"))
writeLines(c(
  "# Validation report: FINN-model-FIA (oregon), zoo v0.2 recipe (stage A + B)", "",
  sprintf("Written %s by `calibration/v02/recipe/oregon/validation.R` (FINN %s, torch %s).", format(Sys.time(), "%Y-%m-%d %H:%M"),
          as.character(packageVersion("FINN")), as.character(packageVersion("torch"))), "",
  "## Sets", "",
  sprintf("- **test**: %d of the %d training plots (`prepared/%s/train_blocked`%s).", length(test_ids), n_train_all, PARAMS$prep_version, if (length(test_ids) < n_train_all) sprintf(", random subset, seed %d", PARAMS$subset_seed) else ""),
  sprintf("- **holdout**: %d plots in 1-degree blocks without any training plot (`prepared/%s/holdout_blocked`); the model has no coordinate predictors, so it is predicted exactly like the test set.",
          nrow(sets$holdout$sites), PARAMS$prep_version), "",
  sprintf("Each plot's 2001-03 stand is simulated over its two %d-year steps (%d test, %d holdout steps); every step is scored at its closing census: states at that census, growth and stem loss as the rate per step, recruitment per ha per step (the aggregation `fit()` uses). Basal area and stems are reported per ha. Persistence (end state = the state at the start of the step) and the change over the step are reported beside the model.",
          PARAMS$step_years, nrow(intervals(sets$test)), nrow(intervals(sets$holdout))), "",
  "Model-side note (audit J): FINN's simulated mortality is floored at 1 % per step (`torch_clamp(m, 0.01, 0.99)` before the death draw); with a decade step this is 0.1 % per year and harmless. No tree table exists for this variant, so the mortality AUC is the site x step x species AUC weighted by the trees at risk. Recruits are drawn with FINN's refined sampler (exact Poisson below a per-patch mean of 2), the only sampler that delivers the small per-patch means of the 0.0168-ha subplots.", "",
  "## Metrics", "", knitr::kable(Mtab), "",
  "## Holdout stand totals by END-census out-of-S basal-area share", "",
  knitr::kable(M_bin[, .(variable, oos_bin_end, n, RMSE = fmt(rmse), R2 = fmt(r2), rmse_persist = fmt(rmse_persist), rmse_ratio_model_persist = fmt(rmse_ratio_model_persist))]), "",
  "## Rates: stem loss and growth per site, calibration ratios", "",
  knitr::kable(AUC[, .(set, site_auc_weighted = fmt(site_auc_weighted), n_rows, n_at_risk, n_died, site_growth_r2 = fmt(site_growth_r2), site_growth_rmse = fmt(site_growth_rmse, 4),
                       deaths_ratio = fmt(deaths_ratio, 2), growth_ratio = fmt(growth_ratio, 2), reg_ratio = fmt(reg_ratio, 2))]), "",
  "The site AUC uses the simulated death probability per step, each site x step x species row contributing `n_died` positives and `n_at_risk - n_died` negatives; the ratios are predicted / observed sums (deaths from the simulated rate, growth per site x species weighted by the trees behind it, recruits per ha per step).", "",
  "## Per species on the holdout (recruitment and death ratios)", "", knitr::kable(SP), "",
  "## Per-species Spearman correlation (holdout)", "",
  knitr::kable(dcast(M_bysp[set == "holdout", .(species_name, variable, spearman = round(spearman, 2))], species_name ~ variable, value.var = "spearman")), "",
  "## Bare-ground succession", "", knitr::kable(EQ[, .(site_lab, ba_last100yr = fmt(ba_last10, 1), ba_prev100yr = fmt(ba_prev10, 1), rel_change = fmt(rel_change), stems_last100yr = fmt(stems_last10, 0))]), "",
  sprintf("%d patches of %.4f ha, %d steps (%d years) from bare ground at six training plots spanning the climate space; `rel_change` compares the mean stand basal area of the last 10 steps (100 years) with the 10 before.", PARAMS$eq_patches, PS, PARAMS$eq_steps, PARAMS$eq_steps * PARAMS$step_years), "",
  "## Figures", "",
  "- `validation/fig1_loss.png` loss per response and stage; `fig2_obs_pred.png` observed vs predicted; `fig3_roc.png` site-level ROC; `fig4_map.png` map of stand basal area; `fig5_bare_ground.png` succession; `fig6_parameters.png` parameter trajectories of both stages.", "",
  "## Fit", "", summ, "",
  "## Provenance", "", knitr::kable(provenance)
), file.path(FIT_DIR, "report.md"))

## ---- 5. hand the card its files -----------------------------------------------------
if (dir.exists(ZOO_DIR)) {
  dir.create(file.path(ZOO_DIR, "figures"), showWarnings = FALSE)
  for (f in c("fig1_loss.png", "fig2_obs_pred.png", "fig3_roc.png", "fig4_map.png", "fig5_bare_ground.png", "fig6_parameters.png"))
    file.copy(file.path(VAL, f), file.path(ZOO_DIR, "figures", f), overwrite = TRUE)
  fwrite(rbind(M[, .(set, level, variable, n, rmse, r2, r, bias, obs_mean)],
               AUC[, .(set, level = "per site x step x species (weighted)", variable = "mortality_auc", n = n_rows, rmse = NA_real_, r2 = NA_real_, r = site_auc_weighted, bias = NA_real_, obs_mean = n_died / n_at_risk)]),
         file.path(ZOO_DIR, "metrics.csv"))
  fwrite(species_dt, file.path(ZOO_DIR, "species.csv"))
}
message("report: ", file.path(FIT_DIR, "report.md"))
