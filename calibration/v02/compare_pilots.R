#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/compare_pilots.R
##
## WHAT   Collects the summaries of every finished v0.2 pilot
##        (calibration/v02/pilots/<pilot>/{summary.csv, metrics.csv}) into one
##        comparison table per variant: holdout R2 of basal area and stems
##        (species x plot and stand total) against persistence, the change
##        metrics, mortality AUC (site level; per tree for conus), the
##        predicted / observed ratios of growth, deaths and recruits, the
##        bare-ground basal area range and drift, runtime.
## RUN    Rscript calibration/v02/compare_pilots.R      (repo root, seconds)
## OUTPUT calibration/v02/PILOTS_v02.md, calibration/v02/pilots_comparison.csv
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))
P <- file.path("calibration", "v02", "pilots")
dirs <- list.dirs(P, recursive = FALSE)
dirs <- dirs[file.exists(file.path(dirs, "summary.csv"))]
if (!length(dirs)) stop("no finished pilot")
rows <- rbindlist(lapply(dirs, function(d) {
  s <- fread(file.path(d, "summary.csv")); m <- fread(file.path(d, "metrics.csv"))
  g <- function(resp, col = "r2") m[response == resp][[col]][1]
  data.table(pilot = basename(d), variant = sub("_.*", "", basename(d)),
             ba_r2 = g("ba species x plot"), ba_r2_persist = g("ba species x plot", "persistence_r2"),
             trees_r2 = g("trees species x plot"), trees_r2_persist = g("trees species x plot", "persistence_r2"),
             ba_tot_r2 = g("ba stand total"), ba_tot_r2_persist = g("ba stand total", "persistence_r2"),
             trees_tot_r2 = g("trees stand total"), trees_tot_r2_persist = g("trees stand total", "persistence_r2"),
             dba_tot_r2 = g("delta ba stand total"), dtrees_tot_r2 = g("delta trees stand total"),
             growth_r2 = g("growth species x plot"), growth_r2_tree = if ("growth_r2_tree" %in% names(s)) s$growth_r2_tree else NA_real_,
             reg_spearman = g("recruitment species x plot", "spearman"),
             mort_auc = s$mort_auc, mort_auc_tree = if ("mort_auc_tree" %in% names(s)) s$mort_auc_tree else NA_real_,
             growth_ratio = s$growth_ratio, deaths_ratio = s$deaths_ratio, reg_ratio = s$reg_ratio, ba_ratio = s$ba_ratio, trees_ratio = s$trees_ratio,
             delta_ba_obs = s$delta_ba_obs, delta_ba_pred = s$delta_ba_pred, delta_trees_obs = s$delta_trees_obs, delta_trees_pred = s$delta_trees_pred,
             bare_ba_min = s$bare_ba_min, bare_ba_max = s$bare_ba_max, bare_drift_max_pct = s$bare_drift_max_pct,
             minutes = s$minutes, epochs = s$epochs, K_max = s$K_max, theta_median = s$theta_median, loss_last = s$loss_last)
}), fill = TRUE)
setorder(rows, variant, pilot)
fwrite(rows, file.path("calibration", "v02", "pilots_comparison.csv"))
num <- names(rows)[vapply(rows, is.numeric, logical(1))]
tab <- copy(rows)[, (num) := lapply(.SD, signif, 3), .SDcols = num]
md <- function(dt) c(paste0("| ", paste(names(dt), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
                     apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
cols1 <- c("pilot", "ba_r2", "ba_r2_persist", "trees_r2", "trees_r2_persist", "ba_tot_r2", "ba_tot_r2_persist", "trees_tot_r2", "trees_tot_r2_persist", "dba_tot_r2", "dtrees_tot_r2")
cols2 <- c("pilot", "growth_r2", "growth_r2_tree", "mort_auc", "mort_auc_tree", "reg_spearman", "growth_ratio", "deaths_ratio", "reg_ratio", "ba_ratio", "trees_ratio")
cols3 <- c("pilot", "delta_ba_obs", "delta_ba_pred", "delta_trees_obs", "delta_trees_pred", "bare_ba_min", "bare_ba_max", "bare_drift_max_pct", "minutes", "epochs", "K_max", "theta_median")
writeLines(c("# v0.2 pilots: holdout comparison", "",
             sprintf("Generated %s by `calibration/v02/compare_pilots.R` from %d pilots. Oregon pilots: full training set (585-590 plots), holdout 118 plots, 300 epochs. CONUS pilots: 300 training plots / 200 holdout plots (fixed seeds), 250 epochs. One thread. `_persist` = the first inventory carried forward; `delta` = change against it; ratios = predicted / observed totals; bare ground = basal area (m2/ha) at the end of a 300-400 year run at three climate points and the largest drift over its last fifth.",
                     format(Sys.time(), "%Y-%m-%d %H:%M"), nrow(rows)), "",
             "## States against persistence", "", md(tab[, ..cols1]), "",
             "## Rates", "", md(tab[, ..cols2]), "",
             "## Change, bare ground, cost", "", md(tab[, ..cols3]), "",
             "## Final losses (baseline-scaled)", "", md(tab[, .(pilot, loss_last)])),
           file.path("calibration", "v02", "PILOTS_v02.md"))
print(tab[, ..cols1]); print(tab[, ..cols2]); print(tab[, ..cols3])
