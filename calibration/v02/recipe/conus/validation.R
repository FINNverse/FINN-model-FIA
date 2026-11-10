#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-zoo | FINN-model-FIA/calibration/v02/recipe/conus/validation.R
##
## WHAT   The validation report of the zoo v0.2 recipe model of the conus variant
##        (fit_stageA.R + fit_stageB.R, fits/conus-annual-v2-recipe/model.pt). Runs
##        on the saved model, so it can be redone after a plot change without a
##        refit. Two sets, named explicitly:
##          TEST    = the plots the model was fitted on (prepared/conus-annual-v2/train,
##                    all 1,500),
##          HOLDOUT = plots the model never saw, in 2-degree blocks that hold no
##                    training plot (prepared/conus-annual-v2/holdout, 1,000).
##        The model has no coordinates as predictors, so the holdout is predicted
##        exactly like the test set.
##        For every plot the first-inventory stand is simulated forward over its
##        own 23 env rows and every observed interval (1-4 per plot, 4-12 years)
##        is aggregated the way fit() scores it: states at the census year, growth
##        and stem loss as the mean annual rate, recruitment as the sum of the
##        annual per-ha rates over the interval. Metrics per response and set
##        (species rows and stand totals) with the PERSISTENCE null (end = the
##        state at the start of the interval: the first inventory, or the previous
##        census) and the change over the interval beside them, change metrics by
##        interval class (audit H), rate calibration ratios (predicted / observed
##        recruits, deaths, growth), holdout errors by the END-census out-of-S
##        share, mortality AUC per tree (FINN::predictTrees) and per site weighted
##        by trees at risk, a per-species Spearman table, and the figures the zoo
##        card shows: loss curves of both stages, observed vs predicted, ROC, map,
##        bare-ground succession across the climate space, parameter trajectories
##        of both stages.
##        Copied from FINN-model-IFN-Spain/calibration/v02/recipe/validation.R
##        (cb41778) and adapted to several intervals per plot, three climate
##        columns and the FIA plot design.
## INPUT  <FIT_DIR>/{model.pt, model_stageA.pt, parameters.pt, loss_history_stageA.csv,
##        loss_history_stageB.csv, fit_summary_stage{A,B}.md}; prepared/conus-annual-v2/{train,holdout}
## OUTPUT <FIT_DIR>/report.md, <FIT_DIR>/validation/*.csv, *.png;
##        copies of metrics.csv, species.csv and the card figures into the zoo
##        folder when it exists (ZOO_DIR below)
## RUN    Rscript calibration/v02/recipe/conus/validation.R        (from the repo root, 4 threads)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
if (dir.exists("~/R/lib-finn-0848b60")) .libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({
  library(FINN); library(data.table); library(torch); library(ggplot2)
})
stopifnot(packageVersion("FINN") == "0.2.0")

DATA_DIR <- file.path("..", "..", "FINN-zoo-data", "fia")
FIT_DIR  <- file.path(DATA_DIR, "fits", "conus-annual-v2-recipe")
VAL      <- file.path(FIT_DIR, "validation")
ZOO_DIR  <- file.path("..", "..", "FINN-zoo-worktrees", "fia", "FIA-conus")   # model-card folder (optional)
dir.create(VAL, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version  = "conus-annual-v2",
  sets          = list(test = "train", holdout = "holdout"),
  n_test_sites  = 1500L,                   # subset of the training plots simulated for the report (all of them)
  subset_seed   = 7L,
  patch_size_ha = 0.01681,                 # as fitted
  n_patches     = 4L,                      # as fitted: the four subplots
  env_vars      = c("tavg_c", "prec_mm", "cwd_z"),   # as fitted (yearly values)
  sim_seed      = 1L,
  interval_classes = c(0, 7, 12, 99),      # change metrics by interval length: <= 7 yr (the eastern 5-7-yr panels), 8-12 yr (the western 10-yr panels), longer
  ## bare-ground succession: plots nearest to the 6 corners / centre of the
  ## (tavg_c, prec_mm) space of the training set (plot means over the env rows),
  ## simulated for eq_years under their mean climate (cwd_z at its plot mean)
  eq_years      = 500L,
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
              init = fread(file.path(d, "init_cohort_dt.csv")), trees = fread(file.path(d, "tree_obs_dt.csv")),
              sites = fread(file.path(d, "siteID_dt.csv")))
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

## ---- 1. simulate every plot over its intervals and aggregate like fit() ----------
## (simulateForest() labels the year driven by env row k as year k; the observation at
## y_end is driven by rows y_end - p + 1 .. y_end)
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
## the state at the START of every interval (persistence null): the first inventory
## for the first interval, the previous census for the later ones (per patch, like obs)
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
## stand totals per interval: ba and trees summed, dbh as the stem-weighted mean; persistence = the start state
TOT <- rbindlist(lapply(names(sets), function(s) {
  tot <- function(x) x[, .(ba = sum(ba), trees = sum(trees), dbh = if (sum(trees) > 0) weighted.mean(dbh, trees, na.rm = TRUE) else NA_real_), by = .(siteID, year)]
  o <- tot(sets[[s]]$obs)[, reg := sets[[s]]$obs[, sum(reg), by = .(siteID, year)]$V1]
  p <- tot(pred[[s]])[, reg := pred[[s]][, sum(reg), by = .(siteID, year)]$V1]
  q <- tot(start[[s]])[, reg := NA_real_]
  m <- merge(melt(o, id.vars = c("siteID", "year"), variable.name = "variable", value.name = "obs"),
             melt(p, id.vars = c("siteID", "year"), variable.name = "variable", value.name = "pred"), by = c("siteID", "year", "variable"))
  m <- merge(m, melt(q, id.vars = c("siteID", "year"), variable.name = "variable", value.name = "persist"), by = c("siteID", "year", "variable"), all.x = TRUE)[, set := s]
  merge(m, intervals(sets[[s]])[, .(siteID, year = y_end, period_length)], by = c("siteID", "year"))
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
M_species <- OP[, metrics(obs, pred), by = .(set, variable)][, level := "species rows (site x interval x species)"]
M_total   <- TOT[, metrics(obs, pred), by = .(set, variable)][, level := "stand totals (site x interval)"]
M_persist <- TOT[variable %in% c("ba", "trees", "dbh"), metrics(obs, persist), by = .(set, variable)][, level := "stand totals: persistence (end = start of the interval)"]
M_change  <- TOT[variable %in% c("ba", "trees"), metrics(obs - persist, pred - persist), by = .(set, variable)][, `:=`(level = "stand totals: change over the interval", variable = paste0("d_", variable))]
M <- rbind(M_species, M_total, M_persist, M_change)
spearman <- function(o, p) { ok <- is.finite(o) & is.finite(p); if (sum(ok) < 3) NA_real_ else suppressWarnings(cor(o[ok], p[ok], method = "spearman")) }
M_bysp <- OP[, .(n = sum(is.finite(obs) & is.finite(pred)), spearman = spearman(obs, pred),
                 bias = mean(pred - obs, na.rm = TRUE), obs_mean = mean(obs, na.rm = TRUE)), by = .(set, variable, species_name)]
fwrite(M, file.path(VAL, "metrics.csv")); fwrite(M_bysp, file.path(VAL, "metrics_by_species.csv"))
## change over the interval by interval class (audit H: the 10-year western panels carry
## the negative observed drift, the 5-7-year eastern ones a positive one)
TOT[, interval_class := cut(period_length, PARAMS$interval_classes, labels = c("<= 7 yr", "8-12 yr", "> 12 yr"))]
M_iv <- TOT[set == "holdout" & variable %in% c("ba", "trees"), c(metrics(obs - persist, pred - persist),
            list(obs_change_mean = mean(obs - persist, na.rm = TRUE), pred_change_mean = mean(pred - persist, na.rm = TRUE))), by = .(variable, interval_class)][order(variable, interval_class)]
fwrite(M_iv, file.path(VAL, "metrics_change_by_interval_class.csv"))
## holdout errors by the END-census out-of-S share (the agreed rare-species rule's stratified evaluation)
BIN <- merge(TOT[set == "holdout" & variable %in% c("ba", "trees")], sets$holdout$sites[, .(siteID, oos_bin_end)], by = "siteID")
M_bin <- BIN[, c(metrics(obs, pred), list(rmse_persist = sqrt(mean((obs - persist)^2)))), by = .(variable, oos_bin_end)][order(variable, oos_bin_end)]
M_bin[, rmse_ratio_model_persist := rmse / rmse_persist]
fwrite(M_bin, file.path(VAL, "metrics_by_oos_bin.csv"))

## ---- 2. rates: AUC per tree and per site, calibration ratios --------------------------
auc_weighted <- function(p, pos, neg) {           # ties count half; pos/neg = counts sharing one probability
  o <- order(p); p <- p[o]; pos <- pos[o]; neg <- neg[o]
  cum_neg <- cumsum(neg) - neg                     # negatives with a smaller p
  ties_neg <- neg                                  # negatives at the same p
  (sum(pos * (cum_neg + 0.5 * ties_neg))) / (sum(pos) * sum(neg))
}
AUC <- rbindlist(lapply(names(sets), function(s) {
  d <- sets[[s]]
  P <- predictTrees(model, tree_data = d$trees, env = d$env, patches = PARAMS$n_patches, patch_size = PS, device = "cpu")
  tr <- copy(d$trees)[, `:=`(growth_pred = P$growth_pred, mort_pred = P$mort_pred)]
  tree_auc <- tr[is.finite(mort_pred), auc_weighted(mort_pred, died, 1 - died)]
  g <- tr[is.finite(growth) & is.finite(growth_pred), metrics(growth, growth_pred)]
  site <- merge(pred[[s]][, .(siteID, year, species, mort)], d$obs[n_at_risk > 0, .(siteID, year, species, n_at_risk, n_died, period_length)],
                by = c("siteID", "year", "species"))
  site[, mort_pred := 1 - (1 - mort)^period_length]            # annual rate -> interval death probability
  site_auc <- site[is.finite(mort_pred), auc_weighted(mort_pred, n_died, n_at_risk - n_died)]
  fwrite(tr, file.path(VAL, sprintf("tree_predictions_%s.csv", s)))
  data.table(set = s, tree_auc = tree_auc, site_auc_weighted = site_auc, n_trees = nrow(tr), n_died = sum(tr$died),
             tree_growth_r2 = g$r2, tree_growth_rmse = g$rmse,
             deaths_ratio = site[is.finite(mort_pred), sum(mort_pred * n_at_risk) / sum(n_died)],
             growth_ratio = tr[is.finite(growth) & is.finite(growth_pred), sum(growth_pred) / sum(growth)],
             reg_ratio = OP[set == s & variable == "reg", sum(pred, na.rm = TRUE) / sum(obs, na.rm = TRUE)])
}))
fwrite(AUC, file.path(VAL, "mortality_auc.csv"))
## per-species recruitment and death ratios on the holdout
SP <- merge(OP[set == "holdout" & variable == "reg", .(reg_obs = sum(obs, na.rm = TRUE), reg_ratio = round(sum(pred, na.rm = TRUE) / sum(obs, na.rm = TRUE), 2)), by = species_name],
            merge(pred$holdout[, .(siteID, year, species, mort)], sets$holdout$obs[n_at_risk > 0, .(siteID, year, species, n_at_risk, n_died, period_length)], by = c("siteID", "year", "species"))[
              , .(n_died = sum(n_died), deaths_ratio = round(sum((1 - (1 - mort)^period_length) * n_at_risk) / sum(n_died), 2)), by = species][
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
OPf <- rbind(OP[, .(set, variable, obs, pred, level = "site x interval x species")], TOT[variable != "reg" | TRUE][, .(set, variable, obs, pred, level = "stand total")])
OPf[, variable := factor(variable, levels = vars, labels = c("basal area (m2/ha)", "stems (/ha)", "mean dbh (cm)", "growth (rel./yr)", "stem loss (/yr)", "recruits (/ha/interval)"))]
g2 <- ggplot(OPf[is.finite(obs) & is.finite(pred)], aes(obs, pred, colour = set)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50") + geom_point(alpha = 0.15, size = 0.6) +
  facet_grid(level ~ variable, scales = "free") + scale_colour_manual(values = c(test = "#1b6ca8", holdout = "#d1495b")) +
  labs(x = "observed", y = "predicted", title = "Observed vs predicted at the census closing each interval", colour = NULL) + theme(legend.position = "bottom")
ggsave(file.path(VAL, "fig2_obs_pred.png"), g2, width = 13, height = 5.5, dpi = 130, bg = "white")
## 3c ROC per tree
roc_dt <- rbindlist(lapply(names(sets), function(s) {
  tr <- fread(file.path(VAL, sprintf("tree_predictions_%s.csv", s)))[is.finite(mort_pred)]
  setorder(tr, -mort_pred)
  tr[, .(tpr = cumsum(died) / sum(died), fpr = cumsum(1 - died) / sum(1 - died), set = s)][seq(1, .N, length.out = min(.N, 2000))]
}))
g3 <- ggplot(roc_dt, aes(fpr, tpr, colour = set)) + geom_abline(slope = 1, intercept = 0, colour = "grey60") + geom_path() +
  scale_colour_manual(values = c(test = "#1b6ca8", holdout = "#d1495b")) + coord_equal() +
  labs(x = "false positive rate", y = "true positive rate", title = "Per-tree stem loss: ROC",
       subtitle = paste(sprintf("%s AUC %.3f", AUC$set, AUC$tree_auc), collapse = "; "), colour = NULL)
ggsave(file.path(VAL, "fig3_roc.png"), g3, width = 5, height = 5, dpi = 130, bg = "white")
## 3d map: observed and predicted stand basal area at the last census of every plot
MAP <- rbindlist(lapply(names(sets), function(s) merge(TOT[set == s & variable == "ba"][, .SD[year == max(year)], by = siteID], sets[[s]]$sites[, .(siteID, lon, lat)], by = "siteID")))
MAPl <- melt(MAP[, .(set, lon, lat, observed = obs, predicted = pred)], id.vars = c("set", "lon", "lat"))
g4 <- ggplot(MAPl, aes(lon, lat, colour = pmin(value, 60))) + geom_point(size = 0.5) + facet_grid(set ~ variable) + coord_quickmap() +
  scale_colour_viridis_c(name = "basal area\n(m2/ha, capped 60)") + labs(title = "Stand basal area at the last census", x = NULL, y = NULL)
ggsave(file.path(VAL, "fig4_map.png"), g4, width = 9, height = 7, dpi = 130, bg = "white")
## 3e bare-ground succession at plots spanning the climate space (plot means of the yearly climate)
train_sites <- merge(sets$test$sites[, .(siteID)], sets$test$env[, lapply(.SD, mean), by = siteID, .SDcols = PARAMS$env_vars], by = "siteID")
pick <- unique(rbindlist(lapply(PARAMS$eq_quantiles, function(q) {
  tq <- quantile(train_sites$tavg_c, q[1]); pq <- quantile(train_sites$prec_mm, q[2])
  z <- train_sites[, (scale(tavg_c)[, 1] - (tq - mean(tavg_c)) / sd(tavg_c))^2 + (scale(prec_mm)[, 1] - (pq - mean(prec_mm)) / sd(prec_mm))^2]
  train_sites[which.min(z)]
})), by = "siteID")
eq_env <- pick[, .(year = seq_len(PARAMS$eq_years)), by = c("siteID", PARAMS$env_vars)][, c("siteID", "year", PARAMS$env_vars), with = FALSE]
eq_env[, siteID := match(siteID, pick$siteID)]
pick[, siteID := seq_len(.N)]
FINN.seed(PARAMS$sim_seed)
eq <- simulateForest(model, env = eq_env, patches = PARAMS$eq_patches, patch_size = PS, init_cohort = NULL, device = "cpu")$long$site
for (cc in c("siteID", "year", "species")) eq[, (cc) := as.integer(as.character(get(cc)))]
eq <- merge(eq[variable %in% c("ba", "trees")], pick, by = "siteID")
eq <- merge(eq, species_dt, by = "species")
eq[, value := value / PS]
eq[, site_lab := sprintf("%.1f C, %.0f mm", tavg_c, prec_mm)]
top <- eq[variable == "ba", .(mx = max(value)), by = .(siteID, species_name)][order(siteID, -mx)][, head(.SD, PARAMS$n_top_species), by = siteID]
eq[, shown := paste(siteID, species_name) %in% paste(top$siteID, top$species_name)]
g5 <- ggplot() +
  geom_line(data = eq[variable == "ba" & !shown], aes(year, value, group = species_name), colour = "grey80", linewidth = 0.3) +
  geom_line(data = eq[variable == "ba" & shown], aes(year, value, colour = species_name), linewidth = 0.6) +
  geom_line(data = eq[variable == "ba", .(value = sum(value)), by = .(site_lab, year)], aes(year, value), colour = "black", linewidth = 0.5, linetype = 2) +
  facet_wrap(~ site_lab, ncol = 3) + labs(x = "year", y = "basal area (m2/ha)", colour = NULL,
  title = "Succession from bare ground (dashed: stand total)", subtitle = sprintf("%d patches of %.4f ha per site under the plot's mean climate; the %d species with the highest basal area coloured", PARAMS$eq_patches, PS, PARAMS$n_top_species)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))
ggsave(file.path(VAL, "fig5_bare_ground.png"), g5, width = 11, height = 7, dpi = 130, bg = "white")
EQ <- eq[variable == "ba", .(ba_total = sum(value)), by = .(siteID, site_lab, year)][
  , .(ba_last100 = mean(ba_total[year > PARAMS$eq_years - 100]), ba_prev100 = mean(ba_total[year > PARAMS$eq_years - 200 & year <= PARAMS$eq_years - 100])), by = .(siteID, site_lab)]
EQ[, rel_change := (ba_last100 - ba_prev100) / ba_prev100]
EQ <- merge(EQ, eq[variable == "trees" & year > PARAMS$eq_years - 100, .(stems_last100 = sum(value) / 100), by = siteID], by = "siteID")
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
  script = c("calibration/conus/01_fiadb_metadata.R .. 07_species_ranking.R", "calibration/v02/species_tradeoff.R",
             "calibration/v02/prepare_conus_v2_annual.R", "calibration/v02/prepare_conus_v2.R",
             "calibration/v02/recipe/conus/fit_stageA.R", "calibration/v02/recipe/conus/fit_stageB.R", "calibration/v02/recipe/conus/validation.R"),
  output = c("FIADB cohort, tree records, WorldClim and TerraClimate per plot, species ranking (calibration/PROVENANCE.md)",
             "calibration/v02/species_set_conus.csv (S = 36, the agreed rare-species rule), species_tradeoff.csv",
             "prepared/annual-v7-S36 (annual step, per-plot intervals, tree table, complete ingrowth target, 2-degree blocks)",
             sprintf("prepared/%s/{train,holdout}, PREPARATION.md, calibration/v02/expected_preparation_conus_v2.csv", PARAMS$prep_version),
             "fits/conus-annual-v2-recipe/{model_stageA.pt, loss_history_stageA.csv, species_parameters_stageA.csv, fit_summary_stageA.md}",
             "fits/conus-annual-v2-recipe/{model.pt, parameters.pt, loss_history_stageB.csv, regeneration_parameters.csv, fit_summary_stageB.md}",
             "fits/conus-annual-v2-recipe/{report.md, validation/}"))
writeLines(c(
  "# Validation report: FINN-model-FIA (conus), zoo v0.2 recipe (stage A + B)", "",
  sprintf("Written %s by `calibration/v02/recipe/conus/validation.R` (FINN %s, torch %s).", format(Sys.time(), "%Y-%m-%d %H:%M"),
          as.character(packageVersion("FINN")), as.character(packageVersion("torch"))), "",
  "## Sets", "",
  sprintf("- **test**: %d of the %d training plots (`prepared/%s/train`%s).", length(test_ids), n_train_all, PARAMS$prep_version, if (length(test_ids) < n_train_all) sprintf(", random subset, seed %d", PARAMS$subset_seed) else ""),
  sprintf("- **holdout**: %d plots in 2-degree blocks without any training plot (`prepared/%s/holdout`); the model has no coordinate predictors, so it is predicted exactly like the test set.",
          nrow(sets$holdout$sites), PARAMS$prep_version), "",
  sprintf("Each plot's first-inventory stand is simulated over its %d env rows; every observed interval (%d test, %d holdout; 1-4 per plot, 4-12 years) is scored at its closing census: states at that year, growth and stem loss as the mean annual rate, recruitment as the sum of the annual per-ha rates (the aggregation `fit()` uses). Basal area and stems are reported per ha. Persistence (end state = the state at the start of the interval) and the change over the interval are reported beside the model.",
          max(sets$test$env$year), nrow(intervals(sets$test)), nrow(intervals(sets$holdout))), "",
  "Model-side note (audit J): FINN's simulated annual mortality is floored at 1 % (`torch_clamp(m, 0.01, 0.99)` before the death draw); the per-tree mortality loss and `predictTrees` see the unfloored rate, the simulated states and the site-level death ratios below the floored one. Recruits are drawn with FINN's refined sampler (exact Poisson below a per-patch mean of 2), the only sampler that delivers the small per-patch means of the 0.0168-ha subplots.", "",
  "## Metrics", "", knitr::kable(Mtab), "",
  "## Holdout change over the interval by interval class", "",
  knitr::kable(M_iv[, .(variable, interval_class, n, RMSE = fmt(rmse), R2 = fmt(r2), r = fmt(r), obs_change_mean = fmt(obs_change_mean), pred_change_mean = fmt(pred_change_mean))]), "",
  "## Holdout stand totals by END-census out-of-S basal-area share", "",
  knitr::kable(M_bin[, .(variable, oos_bin_end, n, RMSE = fmt(rmse), R2 = fmt(r2), rmse_persist = fmt(rmse_persist), rmse_ratio_model_persist = fmt(rmse_ratio_model_persist))]), "",
  "## Rates: stem loss and growth per tree, calibration ratios", "",
  knitr::kable(AUC[, .(set, tree_auc = fmt(tree_auc), site_auc_weighted = fmt(site_auc_weighted), n_trees, n_died, tree_growth_r2 = fmt(tree_growth_r2), tree_growth_rmse = fmt(tree_growth_rmse, 4),
                       deaths_ratio = fmt(deaths_ratio, 2), growth_ratio = fmt(growth_ratio, 2), reg_ratio = fmt(reg_ratio, 2))]), "",
  "Per-tree predictions from `FINN::predictTrees` (the kernels on the observed stand); the site AUC uses the interval death probability implied by the simulated annual rate, each site x interval x species row contributing `n_died` positives and `n_at_risk - n_died` negatives; the ratios are predicted / observed sums (deaths from the simulated rate, growth per tracked tree, recruits per ha per interval).", "",
  "## Per species on the holdout (recruitment and death ratios)", "", knitr::kable(SP), "",
  "## Per-species Spearman correlation (holdout)", "",
  knitr::kable(dcast(M_bysp[set == "holdout", .(species_name, variable, spearman = round(spearman, 2))], species_name ~ variable, value.var = "spearman")), "",
  "## Bare-ground succession", "", knitr::kable(EQ[, .(site_lab, ba_last100 = fmt(ba_last100, 1), ba_prev100 = fmt(ba_prev100, 1), rel_change = fmt(rel_change), stems_last100 = fmt(stems_last100, 0))]), "",
  sprintf("%d patches of %.4f ha, %d years from bare ground at six training plots spanning the climate space, under each plot's mean climate (cwd_z at its plot mean); `rel_change` compares the mean stand basal area of the last 100 years with the 100 before.", PARAMS$eq_patches, PS, PARAMS$eq_years), "",
  "## Figures", "",
  "- `validation/fig1_loss.png` loss per response and stage; `fig2_obs_pred.png` observed vs predicted; `fig3_roc.png` per-tree ROC; `fig4_map.png` map of stand basal area; `fig5_bare_ground.png` succession; `fig6_parameters.png` parameter trajectories of both stages.", "",
  "## Fit", "", summ, "",
  "## Provenance", "", knitr::kable(provenance)
), file.path(FIT_DIR, "report.md"))

## ---- 5. hand the card its files -----------------------------------------------------
if (dir.exists(ZOO_DIR)) {
  dir.create(file.path(ZOO_DIR, "figures"), showWarnings = FALSE)
  for (f in c("fig1_loss.png", "fig2_obs_pred.png", "fig3_roc.png", "fig4_map.png", "fig5_bare_ground.png", "fig6_parameters.png"))
    file.copy(file.path(VAL, f), file.path(ZOO_DIR, "figures", f), overwrite = TRUE)
  fwrite(rbind(M[, .(set, level, variable, n, rmse, r2, r, bias, obs_mean)],
               AUC[, .(set, level = "per tree", variable = "mortality_auc", n = n_trees, rmse = NA_real_, r2 = NA_real_, r = tree_auc, bias = NA_real_, obs_mean = n_died / n_trees)]),
         file.path(ZOO_DIR, "metrics.csv"))
  fwrite(species_dt, file.path(ZOO_DIR, "species.csv"))
}
message("report: ", file.path(FIT_DIR, "report.md"))
