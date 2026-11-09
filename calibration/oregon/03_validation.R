#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/oregon/03_validation.R   (variant "oregon")
##
## WHAT   The validation report of the Oregon model. Rebuilds the model from
##        the released parameter list exactly as the package does
##        (R/definition.R build_model()), so the report validates what a user
##        gets, and scores it on
##          TEST    = the 590 plots it was fitted on (oregon-10yr-v1/train_blocked)
##          HOLDOUT = 118 plots in 1-degree blocks that hold no training plot
##                    (oregon-10yr-v1/holdout_blocked), never seen in fitting.
##        Both are simulated from their own first inventory over the two
##        10-year steps and compared with the two later inventories.
## REPORT per set: pooled metrics for every response (n, RMSE, R2, r, bias,
##        observed mean), species-level and stand-total (summed over species)
##        R2 / RMSE / Spearman for the state variables, mortality AUC weighted
##        by the trees at risk, ROC; figures: loss per response over epochs,
##        observed vs predicted, ROC, map of both sets, succession from bare
##        ground at five points of the climate space (100 steps = 1000 years),
##        parameter trajectories over training; a provenance table.
## RUN    Rscript calibration/oregon/03_validation.R     (from the repository root, ~5 min single-threaded)
## OUTPUT calibration/oregon/report.md, calibration/oregon/metrics.csv,
##        calibration/oregon/validation/*.csv, *.png
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({
  library(FINN); library(data.table); library(torch); library(ggplot2)
})
torch_set_num_threads(1L)

DATA_DIR <- {
  cands <- c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))
  d <- cands[dir.exists(cands)][1]
  if (is.na(d)) stop("FINN-zoo-data/fia not found (expected two levels above the repository, see calibration/README.md)")
  normalizePath(d)
}
VARIANT <- "oregon"
FIT_DIR <- file.path(DATA_DIR, "fits", VARIANT)
CAL_DIR <- file.path("calibration", VARIANT)
VAL     <- file.path(CAL_DIR, "validation"); dir.create(VAL, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  prep_version  = "oregon-10yr-v1",
  test_set      = "train_blocked",
  holdout_set   = "holdout_blocked",
  patch_size_ha = 0.01681,
  n_patches     = 4L,
  step_years    = 10L,                 # one simulated step is one decade
  bare_locations = 5L,                 # points of the climate space for the bare-ground run
  bare_steps    = 100L,                # 100 steps = 1000 years, long enough to settle
  bare_patches  = 15L,
  n_top_species = 8L,                  # coloured in the succession figure; the rest pale grey
  seed          = 1L
)

## ---- helpers -----------------------------------------------------------------
metric <- function(o, p) {
  ok <- is.finite(o) & is.finite(p); o <- o[ok]; p <- p[ok]
  if (length(o) < 3) return(data.table(n = length(o), rmse = NA_real_, r2 = NA_real_, r = NA_real_, spearman = NA_real_, bias = NA_real_, obs_mean = NA_real_))
  sst <- sum((o - mean(o))^2)
  data.table(n = length(o), rmse = sqrt(mean((p - o)^2)),
             r2 = if (sst > 0) 1 - sum((p - o)^2) / sst else NA_real_,
             r = if (sd(o) > 0 && sd(p) > 0) cor(o, p) else NA_real_,
             spearman = if (sd(o) > 0 && sd(p) > 0) cor(o, p, method = "spearman") else NA_real_,
             bias = mean(p) - mean(o), obs_mean = mean(o))
}
## AUC for grouped binary outcomes: each row brings n_pos deaths and n_neg
## survivors sharing one predicted probability; ties count half
group_counts <- function(p, n_pos, n_neg) {
  dt <- data.table(p = p, n_pos = n_pos, n_neg = n_neg)[is.finite(p) & (n_pos + n_neg) > 0]
  dt[, .(n_pos = sum(n_pos), n_neg = sum(n_neg)), by = p][order(p)]
}
weighted_auc <- function(p, n_pos, n_neg) {
  dt <- group_counts(p, n_pos, n_neg)
  if (!nrow(dt) || sum(dt$n_pos) == 0 || sum(dt$n_neg) == 0) return(NA_real_)
  before <- cumsum(dt$n_neg) - dt$n_neg
  sum(dt$n_pos * (before + 0.5 * dt$n_neg)) / (sum(dt$n_pos) * sum(dt$n_neg))
}
roc_points <- function(p, n_pos, n_neg) {
  dt <- group_counts(p, n_pos, n_neg)[order(-p)]
  data.table(fpr = c(0, cumsum(dt$n_neg) / sum(dt$n_neg)), tpr = c(0, cumsum(dt$n_pos) / sum(dt$n_pos)))
}
md_table <- function(dt) {
  dt <- copy(dt); num <- names(dt)[vapply(dt, is.numeric, logical(1))]
  dt[, (num) := lapply(.SD, signif, 3), .SDcols = num]
  c(paste0("| ", paste(names(dt), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
    apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
read_set <- function(set) {
  d <- file.path(DATA_DIR, "prepared", PARAMS$prep_version, set)
  list(obs = fread(file.path(d, "obs_dt.csv")), env = fread(file.path(d, "env_dt.csv")),
       init = fread(file.path(d, "init_cohort_dt.csv")), species = fread(file.path(d, "species_dt.csv"))[order(species)],
       sites = fread(file.path(d, "siteID_dt.csv")))
}

## ---- load: the model as the package rebuilds it ------------------------------
source(file.path("R", "definition.R"))
env_scaling <- fread(file.path(CAL_DIR, "env_scaling.csv"))
model <- build_model(VARIANT, torch_load(file.path(FIT_DIR, "parameters.pt")), env_scaling = as.data.frame(env_scaling))
loss_hist <- fread(file.path(CAL_DIR, "loss_history.csv"))
param_history <- readRDS(file.path(FIT_DIR, "param_history.rds"))
fit_summary <- readLines(file.path(CAL_DIR, "fit_summary.md"))
test <- read_set(PARAMS$test_set); hold <- read_set(PARAMS$holdout_set)
n_species <- nrow(test$species)
message(sprintf("test %d sites | holdout %d sites | %d species", uniqueN(test$obs$siteID), uniqueN(hold$obs$siteID), n_species))

## ---- 1. predictions ----------------------------------------------------------
## One simulation per set from its own first inventory; a simulated step is a
## census, so the states and the per-step rates are compared directly.
predict_set <- function(d) {
  FINN.seed(PARAMS$seed)
  sim <- simulateForest(model, env = d$env, init_cohort = makeInitCohorts(d$init, Nspecies = n_species),
                        patches = PARAMS$n_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")
  S <- dcast(as.data.table(sim$long$site), siteID + year + species ~ variable)
  for (cc in c("siteID", "year", "species")) S[, (cc) := as.integer(as.character(get(cc)))]
  setnames(S, c("dbh", "ba", "trees", "growth", "mort", "r_mean_ha"), c("dbh_m", "ba_m", "trees_m", "growth_m", "mort_m", "reg_m"))
  merge(d$obs, S[, .(siteID, year, species, dbh_m, ba_m, trees_m, growth_m, mort_m, reg_m)], by = c("siteID", "year", "species"))
}
M <- rbind(cbind(set = "test", predict_set(test)), cbind(set = "holdout", predict_set(hold)))
M[, set := factor(set, levels = c("test", "holdout"))]
M[, mort_m := pmin(pmax(mort_m, 0), 1)]
fwrite(M, file.path(VAL, "predictions.csv"))

## ---- 2. metrics --------------------------------------------------------------
vars <- list(
  list("dbh (cm, species mean)",              quote(trees > 0),     "dbh",    "dbh_m"),
  list("basal area (m2 per patch)",           quote(rep(TRUE, .N)), "ba",     "ba_m"),
  list("trees (per patch)",                   quote(rep(TRUE, .N)), "trees",  "trees_m"),
  list("growth (rel. dbh increment per step)", quote(growth_n > 0), "growth", "growth_m"),
  list("mortality (fraction died per step)",  quote(n_at_risk > 0), "mort",   "mort_m"),
  list("recruitment (per ha per step)",       quote(rep(TRUE, .N)), "reg",    "reg_m"))
pooled <- rbindlist(lapply(levels(M$set), function(s) rbindlist(lapply(vars, function(v) {
  d <- M[set == s][eval(v[[2]])]
  cbind(set = s, level = "pooled", variable = v[[1]], species = "all", metric(d[[v[[3]]]], d[[v[[4]]]]))
}))))
by_species <- rbindlist(lapply(levels(M$set), function(s) rbindlist(lapply(vars, function(v) {
  d <- M[set == s][eval(v[[2]])]
  d[, cbind(set = s, level = "species", variable = v[[1]], metric(get(v[[3]]), get(v[[4]]))), by = .(species = species_name)]
}))))
tot <- M[, .(ba = sum(ba), ba_m = sum(ba_m), trees = sum(trees), trees_m = sum(trees_m), reg = sum(reg), reg_m = sum(reg_m)),
         by = .(set, siteID, year)]
total <- rbindlist(lapply(levels(M$set), function(s) rbind(
  cbind(set = s, level = "total", variable = "basal area (m2 per patch)", species = "all (stand total)", metric(tot[set == s]$ba, tot[set == s]$ba_m)),
  cbind(set = s, level = "total", variable = "trees (per patch)",         species = "all (stand total)", metric(tot[set == s]$trees, tot[set == s]$trees_m)),
  cbind(set = s, level = "total", variable = "recruitment (per ha per step)", species = "all (stand total)", metric(tot[set == s]$reg, tot[set == s]$reg_m)))))
auc_tab <- M[n_at_risk > 0, .(n_trees_at_risk = sum(n_at_risk), n_died = sum(n_died),
                              auc = weighted_auc(mort_m, n_died, n_at_risk - n_died)), by = set]
auc_species <- M[n_at_risk > 0, .(n_at_risk = sum(n_at_risk), n_died = sum(n_died),
                                  auc = weighted_auc(mort_m, n_died, n_at_risk - n_died)), by = .(set, species = species_name)]
pooled[variable == "mortality (fraction died per step)", auc := auc_tab$auc[match(set, auc_tab$set)]]
by_species[variable == "mortality (fraction died per step)",
           auc := auc_species$auc[match(paste(set, species), paste(auc_species$set, auc_species$species))]]
metrics <- rbind(pooled, total, by_species, fill = TRUE)[, variant := VARIANT]
setcolorder(metrics, c("variant", "set", "level", "variable", "species"))
fwrite(metrics, file.path(CAL_DIR, "metrics.csv"))
fwrite(auc_tab, file.path(VAL, "mortality_auc.csv"))
print(pooled[, .(set, variable, n, rmse, r2, r, bias, obs_mean, auc)], digits = 3)
print(total[, .(set, variable, n, rmse, r2, r)], digits = 3)

## ---- 3. figures --------------------------------------------------------------
## loss per response
LH <- melt(loss_hist[, .(epoch, dbh, ba, trees, growth, mort, reg)], id.vars = "epoch", variable.name = "response")
g1 <- ggplot(LH, aes(epoch, value, colour = response)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") + geom_line(alpha = 0.85) +
  scale_y_log10() + scale_colour_viridis_d(name = NULL) + theme_minimal(base_size = 10) +
  labs(x = "epoch", y = "loss / intercept-only baseline (1 = no better than the mean)", title = "Training loss per response")
ggsave(file.path(VAL, "fig1_loss_epochs.png"), g1, width = 8, height = 4.5, dpi = 130, bg = "white")

## observed vs predicted
OP <- rbindlist(lapply(vars, function(v) {
  d <- M[eval(v[[2]])]
  data.table(set = d$set, variable = v[[1]], obs = d[[v[[3]]]], pred = d[[v[[4]]]])
}))[is.finite(obs) & is.finite(pred)]
OP[, variable := factor(variable, levels = sapply(vars, `[[`, 1))]
g2 <- ggplot(OP, aes(obs, pred)) +
  geom_hex(bins = 45) + scale_fill_viridis_c(trans = "log10", name = "n") +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = 2, linewidth = 0.4) +
  facet_wrap(~ set + variable, scales = "free", ncol = 6, labeller = labeller(.multi_line = FALSE)) +
  theme_minimal(base_size = 9) + labs(x = "observed", y = "predicted", title = "Observed vs predicted, every response, both sets")
ggsave(file.path(VAL, "fig2_obs_pred_scatter.png"), g2, width = 15, height = 6, dpi = 130, bg = "white")

## ROC
RC <- M[n_at_risk > 0, roc_points(mort_m, n_died, n_at_risk - n_died), by = set]
RC <- merge(RC, auc_tab[, .(set, label = sprintf("%s  (AUC %.3f)", set, auc))], by = "set")
g3 <- ggplot(RC, aes(fpr, tpr, colour = label)) + geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
  geom_path(linewidth = 0.7) + coord_equal() + theme_minimal(base_size = 10) +
  labs(x = "false positive rate", y = "true positive rate", colour = NULL, title = "Mortality ROC, trees at risk weighted") +
  theme(legend.position = "bottom")
ggsave(file.path(VAL, "fig3_mortality_roc.png"), g3, width = 5.5, height = 6, dpi = 130, bg = "white")

## map: observed and predicted stand basal area at the last census
last <- tot[, .SD[year == max(year)], by = .(set, siteID)]
MP <- rbind(merge(last[set == "test"], test$sites[, .(siteID, lon, lat)], by = "siteID"),
            merge(last[set == "holdout"], hold$sites[, .(siteID, lon, lat)], by = "siteID"))
MP <- melt(MP[, .(set, lon, lat, observed = ba / PARAMS$patch_size_ha, predicted = ba_m / PARAMS$patch_size_ha)],
           id.vars = c("set", "lon", "lat"), variable.name = "what", value.name = "ba_ha")
g4 <- ggplot(MP, aes(lon, lat, colour = ba_ha, shape = set)) + geom_point(size = 1.6, alpha = 0.85) +
  scale_colour_viridis_c(name = expression("basal area m"^2*" ha"^-1), limits = c(0, quantile(MP$ba_ha, 0.98)), oob = scales::squish) +
  scale_shape_manual(values = c(test = 16, holdout = 17)) + facet_wrap(~ what) + coord_quickmap() + theme_minimal(base_size = 10) +
  labs(x = NULL, y = NULL, title = "Stand basal area at the last inventory: training plots (circles) and holdout blocks (triangles)")
ggsave(file.path(VAL, "fig4_map.png"), g4, width = 12, height = 5.5, dpi = 130, bg = "white")

## succession from bare ground at five points of the climate space
site_env <- test$env[, .(temp = mean(temp), prec = mean(prec)), by = siteID][order(temp)]
picks <- site_env[round(quantile(seq_len(.N), c(0.1, 0.3, 0.5, 0.7, 0.9)))]
picks[, loc := factor(sprintf("%.1f C, %.0f mm", temp, prec), levels = sprintf("%.1f C, %.0f mm", temp, prec))][, sim_site := .I]
env_bare <- picks[, .(year = seq_len(PARAMS$bare_steps)), by = .(siteID = sim_site, temp, prec)]
setorder(env_bare, siteID, year)
FINN.seed(PARAMS$seed)
sim_bare <- simulateForest(model, env = env_bare[, .(siteID, year, temp, prec)], init_cohort = NULL,
                           patches = PARAMS$bare_patches, patch_size = PARAMS$patch_size_ha, device = "cpu")
SB <- as.data.table(sim_bare$long$site)[variable %in% c("ba", "trees")]
for (cc in c("siteID", "year", "species")) SB[, (cc) := as.integer(as.character(get(cc)))]
SB <- dcast(SB, siteID + year + species ~ variable)
SB[, `:=`(ba_ha = ba / PARAMS$patch_size_ha, trees_ha = trees / PARAMS$patch_size_ha, years = year * PARAMS$step_years)]
SB <- merge(SB, test$species[, .(species, species_name)], by = "species")
SB <- merge(SB, picks[, .(siteID = sim_site, loc)], by = "siteID")
fwrite(SB[, .(loc, step = year, years, species_name, ba_ha, trees_ha)], file.path(VAL, "bare_ground_succession.csv"))
top_species <- SB[, .(peak = max(ba_ha)), by = species_name][order(-peak)][seq_len(PARAMS$n_top_species), species_name]
g5 <- ggplot() +
  geom_line(data = SB[!species_name %in% top_species], aes(years, ba_ha, group = species_name), colour = "grey82", linewidth = 0.3) +
  geom_line(data = SB[species_name %in% top_species], aes(years, ba_ha, colour = species_name), linewidth = 0.6) +
  facet_wrap(~ loc, ncol = 5) + scale_colour_brewer(palette = "Paired", name = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom") +
  labs(x = "years since bare ground", y = expression("basal area  m"^2*" ha"^-1),
       title = sprintf("Succession from bare ground at five points of the climate space (%d steps of %d years, %d patches)",
                       PARAMS$bare_steps, PARAMS$step_years, PARAMS$bare_patches),
       subtitle = sprintf("the %d species reaching the highest basal area at any time are coloured; every other species is one pale grey line", PARAMS$n_top_species))
ggsave(file.path(VAL, "fig5_bare_ground_succession.png"), g5, width = 15, height = 5.5, dpi = 130, bg = "white")
eq_window <- PARAMS$bare_steps * c(0.8, 0.9, 1)
equil <- SB[, .(ba_last = sum(ba_ha[year > eq_window[2]]) / sum(year > eq_window[2]) * n_species,
                ba_prev = sum(ba_ha[year > eq_window[1] & year <= eq_window[2]]) / sum(year > eq_window[1] & year <= eq_window[2]) * n_species,
                trees_last = sum(trees_ha[year > eq_window[2]]) / sum(year > eq_window[2]) * n_species,
                n_species_above_1 = uniqueN(species_name[year > eq_window[2] & ba_ha > 1])), by = loc]
equil[, drift_pct := 100 * (ba_last - ba_prev) / ba_prev]
fwrite(equil, file.path(VAL, "bare_ground_equilibrium.csv"))
print(equil, digits = 3)

## parameter trajectories
if (length(param_history) > 1) {
  step <- nrow(loss_hist) / length(param_history)
  PT <- rbindlist(lapply(seq_along(param_history), function(i) {
    s <- param_history[[i]]
    rbindlist(lapply(names(s), function(nm) {
      mat <- as.matrix(s[[nm]])
      data.table(epoch = i * step, block = nm, row = as.integer(row(mat)), col = as.integer(col(mat)), value = as.numeric(mat))
    }))
  }))
  PT[, kind := fifelse(grepl("^nn_", block), "climate response", "mechanistic")]
  fwrite(PT, file.path(VAL, "parameter_trajectories.csv"))
  mech <- PT[kind == "mechanistic"]
  mech[, panel := fifelse(block == "reg_logK", "reg_logK -> K (recruits/ha/step)", sprintf("%s [, %d]", sub("_unconstrained", "", block), col))]
  mech[block == "reg_logK", value := exp(value)]
  g6a <- ggplot(mech, aes(epoch, value, group = row)) + geom_line(alpha = 0.35, linewidth = 0.3, colour = "#1b6ca8") +
    facet_wrap(~ panel, scales = "free_y", ncol = 4) + theme_minimal(base_size = 9) +
    labs(x = "epoch", y = "value (unconstrained scale)", title = "Mechanistic per-species parameters over training",
         subtitle = sprintf("one line per species; %d snapshots, every %g epochs", length(param_history), step))
  setorder(PT, block, row, col, epoch)
  PT[, delta := abs(value - shift(value)), by = .(block, row, col)]
  conv <- PT[!is.na(delta), .(mean_abs_change = mean(delta)), by = .(epoch, block, kind)]
  g6b <- ggplot(conv, aes(epoch, mean_abs_change, colour = block)) + geom_line(linewidth = 0.5) + geom_point(size = 1) + scale_y_log10() +
    facet_wrap(~ kind, ncol = 2, scales = "free_y") + theme_minimal(base_size = 9) +
    theme(legend.position = "right", legend.text = element_text(size = 6), legend.key.height = unit(7, "pt")) +
    labs(x = "epoch", y = "mean |change| since the previous snapshot (log)", title = "Are the parameters still moving?",
         subtitle = "falling towards zero means the fit settled; a flat line means it did not")
  ggsave(file.path(VAL, "fig6_parameters.png"), gridExtra::grid.arrange(g6a, g6b, ncol = 1, heights = c(1.1, 1)),
         width = 12, height = 10, dpi = 130, bg = "white")
}

## ---- 4. report.md ------------------------------------------------------------
sep <- fread(file.path(DATA_DIR, "prepared", PARAMS$prep_version, "holdout_separation.csv"))
sp_state <- dcast(by_species[variable %in% c("basal area (m2 per patch)", "trees (per patch)") & set == "holdout",
                             .(species, variable = fifelse(grepl("basal", variable), "ba", "trees"), r2, rmse, spearman)],
                  species ~ variable, value.var = c("r2", "rmse", "spearman"))
sp_state <- merge(sp_state, by_species[variable == "mortality (fraction died per step)" & set == "holdout", .(species, mort_auc = auc)], by = "species", all.x = TRUE)
sp_state <- merge(sp_state, by_species[variable == "growth (rel. dbh increment per step)" & set == "holdout", .(species, growth_r2 = r2)], by = "species", all.x = TRUE)
sp_state <- merge(sp_state, by_species[variable == "recruitment (per ha per step)" & set == "holdout", .(species, reg_spearman = spearman)], by = "species", all.x = TRUE)
verdict_ok <- all(is.finite(equil$ba_last)) && all(equil$ba_last > 5 & equil$ba_last < 150) && all(abs(equil$drift_pct) < 10)
writeLines(c(
  sprintf("# Validation report: FINN-model-FIA variant `%s`", VARIANT), "",
  sprintf("Generated %s by `calibration/%s/03_validation.R` (FINN %s, torch %s). The model is rebuilt from `parameters.pt` with `build_model()`, as the package does.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), VARIANT, as.character(packageVersion("FINN")), as.character(packageVersion("torch"))), "",
  fit_summary[-1], "",
  "## The two sets", "",
  sprintf("- **test**: the %d plots the model was fitted on (`%s/%s`).", uniqueN(test$obs$siteID), PARAMS$prep_version, PARAMS$test_set),
  sprintf("- **holdout**: %d plots in 1-degree blocks that hold no training plot (`%s/%s`), never used in fitting. Distance to the nearest training plot: median %.0f km, minimum %.0f km, 10 %% below %.0f km.",
          uniqueN(hold$obs$siteID), PARAMS$prep_version, PARAMS$holdout_set, median(sep$nearest_train_km), min(sep$nearest_train_km), quantile(sep$nearest_train_km, 0.1)), "",
  "Each set is simulated from its own first inventory over two 10-year steps; the state at a step is compared with the inventory of that step, the per-step rates with the observed rates of the interval that ends there.", "",
  "## Pooled metrics (all species x plots x steps)", "",
  md_table(pooled[, .(set, variable, n, rmse, r2, r, bias, obs_mean, auc)]), "",
  "## Stand totals (summed over species per plot and step)", "",
  md_table(total[, .(set, variable, n, rmse, r2, r, bias, obs_mean)]), "",
  "## Mortality AUC (trees at risk weighted)", "", md_table(auc_tab), "",
  "## Holdout, per species", "",
  "R2, RMSE and Spearman of basal area and stem number per patch, the per-tree mortality AUC, growth R2 and recruitment Spearman, each within the species.", "",
  md_table(sp_state[order(-r2_ba)]), "",
  "## Succession from bare ground", "",
  sprintf("Five climate points (10/30/50/70/90 %% quantile of plot temperature), %d steps of %d years, %d patches. Stand basal area (m2/ha) averaged over the last 10 %% of the run against the 10 %% before, and the number of species above 1 m2/ha.",
          PARAMS$bare_steps, PARAMS$step_years, PARAMS$bare_patches), "",
  md_table(equil[, .(loc, ba_last, ba_prev, drift_pct, trees_last, n_species_above_1)]), "",
  sprintf("Verdict: %s", if (verdict_ok) "the stand settles at a plausible basal area (5-150 m2/ha) with less than 10 % drift over the last fifth of the run at every climate point." else "CHECK: at least one climate point drifts, collapses or runs away (see the table); the model is not validated as an equilibrium model there."), "",
  "## Figures", "",
  "![loss](validation/fig1_loss_epochs.png)", "![observed vs predicted](validation/fig2_obs_pred_scatter.png)",
  "![ROC](validation/fig3_mortality_roc.png)", "![map](validation/fig4_map.png)",
  "![bare ground](validation/fig5_bare_ground_succession.png)", "![parameters](validation/fig6_parameters.png)", "",
  "## Provenance", "",
  "| step | script | inputs | outputs |", "|---|---|---|---|",
  "| FIA raw | USDA FS FIA DataMart, CSV dump (June 2025) | - | ENTIRE_TREE/PLOT/COND/SUBPLOT csv |",
  "| Oregon export | FINN-fia scripts/01-03, 07 (2026-06; FINN commit 15922a1 data-raw) | FIADB csv, WorldClim 2.1 | full_tree_dt.csv, env_unscaled_dt.csv, species_dt.csv |",
  "| coordinates | calibration/oregon/01_prepare.R | FIADB ENTIRE_PLOT.csv | oregon_plot_coords.csv |",
  sprintf("| preparation | calibration/oregon/01_prepare.R | the export + coordinates | prepared/%s/{all,train_blocked,holdout_blocked}, PREPARATION.md, expected_preparation.csv |", PARAMS$prep_version),
  sprintf("| fit | calibration/oregon/02_fit.R | prepared/%s/%s | fits/oregon/parameters.pt, model.pt, loss_history.csv, fit_summary.md, env_scaling.csv |", PARAMS$prep_version, PARAMS$test_set),
  sprintf("| validation | calibration/oregon/03_validation.R | parameters.pt, prepared/%s/{%s,%s} | this report, metrics.csv, validation/ |", PARAMS$prep_version, PARAMS$test_set, PARAMS$holdout_set)
), file.path(CAL_DIR, "report.md"))
message("report: ", file.path(CAL_DIR, "report.md"))
