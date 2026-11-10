## FINN-model-FIA | calibration/v02/pilot_eval.R
##
## Generic evaluation helpers sourced by every v0.2 pilot script (metric
## functions only; nothing here changes what a pilot fits). pilot_evaluate()
## simulates a fitted model on a holdout set, aggregates the annual predictions
## over the observed intervals as fit() scores them, and writes one metrics
## table with, per response, the model against the persistence baseline (the
## first inventory carried forward), the change metrics (delta basal area, delta
## stems against the first inventory, species x plot and stand total), the
## mortality AUC weighted by trees at risk, the recruitment ratio, and a short
## bare-ground run (basal area at the end, drift over the last fifth).
suppressPackageStartupMessages({ library(FINN); library(data.table); library(torch) })

pe_metric <- function(o, p) {
  ok <- is.finite(o) & is.finite(p); o <- o[ok]; p <- p[ok]
  if (length(o) < 3) return(data.table(n = length(o), rmse = NA_real_, r2 = NA_real_, spearman = NA_real_, bias = NA_real_))
  sst <- sum((o - mean(o))^2)
  data.table(n = length(o), rmse = sqrt(mean((p - o)^2)), r2 = if (sst > 0) 1 - sum((p - o)^2) / sst else NA_real_,
             spearman = if (sd(o) > 0 && sd(p) > 0) cor(o, p, method = "spearman") else NA_real_, bias = mean(p) - mean(o))
}
pe_weighted_auc <- function(p, n_pos, n_neg) {
  dt <- data.table(p = p, n_pos = n_pos, n_neg = n_neg)[is.finite(p) & (n_pos + n_neg) > 0]
  dt <- dt[, .(n_pos = sum(n_pos), n_neg = sum(n_neg)), by = p][order(p)]
  if (!nrow(dt) || sum(dt$n_pos) == 0 || sum(dt$n_neg) == 0) return(NA_real_)
  before <- cumsum(dt$n_neg) - dt$n_neg
  sum(dt$n_pos * (before + 0.5 * dt$n_neg)) / (sum(dt$n_pos) * sum(dt$n_neg))
}
pe_read_set <- function(d) {
  list(obs = fread(file.path(d, "obs_dt.csv")), env = fread(file.path(d, "env_dt.csv")),
       init = fread(file.path(d, "init_cohort_dt.csv")), species = fread(file.path(d, "species_dt.csv"))[order(species)],
       trees = if (file.exists(file.path(d, "tree_obs_dt.csv"))) fread(file.path(d, "tree_obs_dt.csv")) else NULL)
}

## model: a fitted finn model; d: a set as read by pe_read_set(); step_years: 10 (oregon) or 1 (conus)
## env_vars: the covariate names; bare_env: a data.table of a few env rows (one per location) for the bare-ground run
pilot_evaluate <- function(model, d, label, out_dir, n_patches = 4L, patch_size = 0.01681, step_years = 1L,
                           bare_env = NULL, bare_years = 300L, bare_patches = 10L, seed = 1L, batchsize = 250L) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  model$eval()
  n_species <- nrow(d$species)
  init <- d$init
  if ("living" %in% names(init)) init <- init[living == TRUE]
  FINN.seed(seed)
  sim <- simulateForest(model, env = d$env, init_cohort = makeInitCohorts(d$init, Nspecies = n_species),
                        patches = n_patches, patch_size = patch_size, batchsize = batchsize, device = "cpu")
  S <- dcast(as.data.table(sim$long$site), siteID + year + species ~ variable)
  for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
  obs <- copy(d$obs)
  if (!"period_length" %in% names(obs)) obs[, period_length := 1L]
  O <- obs[, .(siteID, species, y_end = year, y0 = year - period_length)]; O[, rid := .I]
  W <- S[O, on = .(siteID, species, year > y0, year <= y_end), nomatch = 0L, allow.cartesian = TRUE,
         .(rid = i.rid, growth = x.growth, mort = x.mort, r_mean_ha = x.r_mean_ha)
        ][, .(growth_m = mean(growth, na.rm = TRUE), mort_m = mean(mort, na.rm = TRUE), reg_m = sum(r_mean_ha)), by = rid]
  St <- merge(O, S[, .(siteID, species, y_end = year, dbh_m = dbh, ba_m = ba, trees_m = trees)], by = c("siteID", "species", "y_end"), all.x = TRUE)
  P <- merge(St, W, by = "rid", all.x = TRUE)[, .(siteID, species, year = y_end, dbh_m, ba_m, trees_m, growth_m, mort_m, reg_m)]
  M <- merge(obs, P, by = c("siteID", "year", "species"))
  M[, mort_m := pmin(pmax(mort_m, 0), 1)]
  M[, mort_period_m := if (step_years > 1) mort_m else 1 - (1 - mort_m)^period_length]
  i0 <- init[trees > 0, .(ba_p = sum(pi * (dbh / 200)^2 * trees) / n_patches, trees_p = sum(trees) / n_patches), by = .(siteID, species)]
  M <- merge(M, i0, by = c("siteID", "species"), all.x = TRUE)
  M[is.na(ba_p), ba_p := 0]; M[is.na(trees_p), trees_p := 0]
  T <- M[, .(ba = sum(ba), ba_m = sum(ba_m), ba_p = sum(ba_p), trees = sum(trees), trees_m = sum(trees_m), trees_p = sum(trees_p)), by = .(siteID, year)]
  rows <- list(
    cbind(response = "ba species x plot", model = pe_metric(M$ba, M$ba_m), persistence_r2 = pe_metric(M$ba, M$ba_p)$r2),
    cbind(response = "trees species x plot", model = pe_metric(M$trees, M$trees_m), persistence_r2 = pe_metric(M$trees, M$trees_p)$r2),
    cbind(response = "ba stand total", model = pe_metric(T$ba, T$ba_m), persistence_r2 = pe_metric(T$ba, T$ba_p)$r2),
    cbind(response = "trees stand total", model = pe_metric(T$trees, T$trees_m), persistence_r2 = pe_metric(T$trees, T$trees_p)$r2),
    cbind(response = "delta ba species x plot", model = pe_metric(M$ba - M$ba_p, M$ba_m - M$ba_p), persistence_r2 = NA_real_),
    cbind(response = "delta ba stand total", model = pe_metric(T$ba - T$ba_p, T$ba_m - T$ba_p), persistence_r2 = NA_real_),
    cbind(response = "delta trees stand total", model = pe_metric(T$trees - T$trees_p, T$trees_m - T$trees_p), persistence_r2 = NA_real_),
    cbind(response = "growth species x plot", model = pe_metric(M[growth_n > 0]$growth, M[growth_n > 0]$growth_m), persistence_r2 = NA_real_),
    cbind(response = "mortality species x plot", model = pe_metric(M[n_at_risk > 0]$mort, M[n_at_risk > 0]$mort_m), persistence_r2 = NA_real_),
    cbind(response = "recruitment species x plot", model = pe_metric(M$reg, M$reg_m), persistence_r2 = NA_real_))
  out <- rbindlist(rows)
  setnames(out, sub("^model\\.", "", names(out)))
  out[, pilot := label]
  extra <- data.table(pilot = label,
    mort_auc = M[n_at_risk > 0, pe_weighted_auc(mort_period_m, n_died, n_at_risk - n_died)],
    growth_ratio = M[growth_n > 0, weighted.mean(growth_m, growth_n, na.rm = TRUE) / weighted.mean(growth, growth_n, na.rm = TRUE)],
    deaths_ratio = M[, sum(mort_period_m * n_at_risk, na.rm = TRUE) / sum(n_died)],
    reg_ratio = M[, sum(reg_m) / sum(reg)], ba_ratio = M[, sum(ba_m) / sum(ba)], trees_ratio = M[, sum(trees_m) / sum(trees)],
    delta_ba_obs = mean(T$ba - T$ba_p), delta_ba_pred = mean(T$ba_m - T$ba_p),
    delta_trees_obs = mean(T$trees - T$trees_p), delta_trees_pred = mean(T$trees_m - T$trees_p))
  ## per-tree AUC when a tree table is there (conus)
  if (!is.null(d$trees)) {
    TRp <- predictTrees(model, tree_data = d$trees, env = d$env, patches = n_patches, patch_size = patch_size, device = "cpu")
    extra[, mort_auc_tree := pe_weighted_auc(TRp$mort_pred, d$trees$died, 1 - d$trees$died)]
    ok <- is.finite(d$trees$growth) & is.finite(TRp$growth_pred)
    extra[, growth_r2_tree := pe_metric(d$trees$growth[ok], TRp$growth_pred[ok])$r2]
  }
  ## bare ground
  if (!is.null(bare_env)) {
    eb <- bare_env[, .(year = seq_len(bare_years)), by = c("siteID", setdiff(names(bare_env), c("siteID", "year")))]
    setcolorder(eb, c("siteID", "year")); setorder(eb, siteID, year)
    FINN.seed(seed)
    sb <- simulateForest(model, env = eb, init_cohort = NULL, patches = bare_patches, patch_size = patch_size, device = "cpu")
    B <- as.data.table(sb$long$site)[variable == "ba"]
    for (cc in c("siteID", "year")) B[, (cc) := as.integer(as.character(get(cc)))]
    Bt <- B[, .(ba_ha = sum(value) / patch_size), by = .(siteID, year)]
    eq <- Bt[, .(ba_last = mean(ba_ha[year > 0.9 * bare_years]), ba_prev = mean(ba_ha[year > 0.8 * bare_years & year <= 0.9 * bare_years])), by = siteID]
    extra[, `:=`(bare_ba_min = min(eq$ba_last), bare_ba_max = max(eq$ba_last), bare_drift_max_pct = max(abs(100 * (eq$ba_last - eq$ba_prev) / eq$ba_prev)))]
    fwrite(Bt, file.path(out_dir, "bare_ground.csv"))
  }
  fwrite(out, file.path(out_dir, "metrics.csv")); fwrite(extra, file.path(out_dir, "summary.csv")); fwrite(M, file.path(out_dir, "predictions.csv"))
  print(out[, .(response, n, rmse = signif(rmse, 3), r2 = signif(r2, 3), spearman = signif(spearman, 3), persistence_r2 = signif(persistence_r2, 3))])
  print(extra)
  invisible(list(metrics = out, summary = extra, predictions = M))
}

## observed per-species recruitment cap: the q-quantile of the observed recruits per
## ha per period on the training set, times a factor, with a floor (per ha per STEP)
pe_reg_cap <- function(obs, species, q = 0.95, factor = 2, floor_cap = 20, per_step = 1) {
  caps <- obs[, .(cap = pmax(floor_cap, factor * quantile(reg, q, na.rm = TRUE)) / per_step), by = species]
  caps <- merge(data.table(species = species$species), caps, by = "species", all.x = TRUE)
  caps[is.na(cap), cap := floor_cap / per_step]
  caps[order(species)]$cap
}
