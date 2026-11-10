#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/audit_v02.R
##
## WHAT   The v0.2 audit of the data and of the v0.1 fits of both variants, the
##        numbers behind calibration/AUDIT_v02.md. For each variant:
##          A  env rows vs observation windows (trap 1): steps before the first
##             census == first period_length; env covers every census
##          B  ingrowth vs recruitment (trap 2): how the FIA design tallies a
##             tree crossing 12.7 cm (microplot record before -> PREV_TRE_CN set
##             -> status "alive", NOT "new"), so what share of the true ingrowth
##             the v0.1 recruitment target missed; RECONCILECD of the "new" trees
##          C  harvest and the mortality target: cut trees and land-use deaths
##             left in the cohort (M0/M1 tiers)
##          D  macroplot thinning: how many trees and which plots it touches
##          E  residuals of the v0.1 holdout predictions: rate calibration
##             (predicted / observed growth, mortality, recruitment) per species,
##             change metrics (delta BA, delta stems against the first inventory)
##             per species and stand, by interval length and by initial size
##          F  convergence: mean absolute parameter change per block over the
##             last snapshots of param_history; recruitment dispersion theta; K
##          G  the "other" class (oregon) in the fit and in the bare-ground run
## INPUT  <DATA_DIR>/prepared/{oregon-10yr-v1,conus-annual-v1}, intermediate/conus
##        (calibration_set_trees.csv), fits/<variant>/{param_history.rds,parameters.pt},
##        calibration/<variant>/{validation/predictions.csv, validation/tree_predictions.csv,
##        validation/bare_ground_succession.csv, loss_history.csv}
## OUTPUT calibration/v02/audit/*.csv, *.png
## RUN    Rscript calibration/v02/audit_v02.R     (from the repository root, ~5 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(torch); library(FINN); library(ggplot2) })
torch_set_num_threads(1L)
DATA_DIR <- {
  cands <- c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))
  d <- cands[dir.exists(cands)][1]
  if (is.na(d)) stop("FINN-zoo-data/fia not found")
  normalizePath(d)
}
OUT <- file.path("calibration", "v02", "audit"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PARAMS <- list(min_dbh_cm = 12.7, macro_dbh = 61, patch_size_ha = 0.01681, n_patches = 4L)
metric <- function(o, p) {
  ok <- is.finite(o) & is.finite(p); o <- o[ok]; p <- p[ok]
  if (length(o) < 3) return(data.table(n = length(o), rmse = NA_real_, r2 = NA_real_, spearman = NA_real_, bias = NA_real_))
  sst <- sum((o - mean(o))^2)
  data.table(n = length(o), rmse = sqrt(mean((p - o)^2)), r2 = if (sst > 0) 1 - sum((p - o)^2) / sst else NA_real_,
             spearman = if (sd(o) > 0 && sd(p) > 0) cor(o, p, method = "spearman") else NA_real_, bias = mean(p) - mean(o))
}
sets <- list(oregon = list(prep = "oregon-10yr-v1", train = "train_blocked", holdout = "holdout_blocked", step = 10),
             conus  = list(prep = "conus-annual-v1", train = "train", holdout = "holdout", step = 1))
lines <- character(0); say <- function(...) { s <- sprintf(...); message(s); lines <<- c(lines, s) }

## ---- A. alignment ------------------------------------------------------------
for (v in names(sets)) {
  s <- sets[[v]]
  for (st in c("train", "holdout")) {
    d <- file.path(DATA_DIR, "prepared", s$prep, s[[st]])
    obs <- fread(file.path(d, "obs_dt.csv")); env <- fread(file.path(d, "env_dt.csv"))
    first <- obs[, .(first_year = min(year), p1 = if ("period_length" %in% names(obs)) period_length[which.min(year)] else 1L), by = siteID]
    envc <- env[, .(env_min = min(year), env_max = max(year), n_env = .N), by = siteID]
    chk <- merge(first, envc, by = "siteID")
    last <- obs[, .(last = max(year)), by = siteID]; chk <- merge(chk, last, by = "siteID")
    say("[A %s/%s] sites %d; env starts at 1: %d; env rows == max env year: %d; first census year == first period_length: %d; env covers last census: %d",
        v, st, nrow(chk), sum(chk$env_min == 1), sum(chk$n_env == chk$env_max), sum(chk$first_year == chk$p1), sum(chk$env_max >= chk$last))
  }
}

## ---- B. ingrowth vs recruitment (conus cohort trees) ------------------------
tr <- fread(file.path(DATA_DIR, "intermediate", "conus", "calibration_set_trees.csv"),
            select = c("plotID", "SUBP", "TREE", "CN", "PREV_TRE_CN", "INVYR", "SPCD", "STATUSCD", "RECONCILECD", "geom", "dbh_cm", "TPA_UNADJ"))
sites <- unique(rbind(fread(file.path(DATA_DIR, "prepared", "conus-annual-v1", "train", "siteID_dt.csv"))[, .(siteName)],
                      fread(file.path(DATA_DIR, "prepared", "conus-annual-v1", "holdout", "siteID_dt.csv"))[, .(siteName)]))
tr <- tr[plotID %chin% sites$siteName]
tr[, treeName := paste(plotID, SUBP, TREE, sep = "_")]
setorder(tr, treeName, INVYR)
tr[, `:=`(dbh_prev = shift(dbh_cm), geom_prev = shift(geom), status_prev = shift(STATUSCD), invyr_prev = shift(INVYR)), by = treeName]
tr[, first_plot_year := min(INVYR), by = plotID]
cand <- tr[STATUSCD == 1 & dbh_cm >= PARAMS$min_dbh_cm & INVYR > first_plot_year & geom == "subp"]
cand[, kind := fcase(is.na(PREV_TRE_CN) | PREV_TRE_CN == "", "no previous CN (counted as recruit in v0.1)",
                     !is.na(dbh_prev) & dbh_prev < PARAMS$min_dbh_cm, "tracked from below 12.7 cm (microplot record; NOT counted in v0.1)",
                     !is.na(dbh_prev) & dbh_prev >= PARAMS$min_dbh_cm, "survivor >= 12.7 cm before",
                     default = "previous CN but no previous record in the extract")]
B1 <- cand[, .N, by = kind][order(-N)]
B1[, share := N / sum(N)]
fwrite(B1, file.path(OUT, "B_ingrowth_kinds_conus.csv"))
say("[B conus] trees >= 12.7 cm alive at a later census on the 2500 model plots: %s", paste(sprintf("%s: %d", B1$kind, B1$N), collapse = "; "))
ing <- cand[kind != "survivor >= 12.7 cm before"]
B2 <- ing[, .N, by = .(kind, RECONCILECD)][order(kind, -N)]
fwrite(B2, file.path(OUT, "B_ingrowth_reconcilecd_conus.csv"))
n_new <- ing[grepl("no previous", kind), .N]; n_tracked <- ing[grepl("tracked", kind), .N]
say("[B conus] v0.1 recruitment target counts %d trees; %d further trees crossed 12.7 cm with a microplot record before (missed): target is %.0f %% of the ingrowth", n_new, n_tracked, 100 * n_new / (n_new + n_tracked))
## microplot share of the trees 2.5-12.7 cm (what feeds the ingrowth) and the ratio to the subplot area
B3 <- tr[STATUSCD == 1, .(n = .N), by = .(geom, below = dbh_cm < PARAMS$min_dbh_cm)][order(geom, below)]
fwrite(B3, file.path(OUT, "B_geom_by_size_conus.csv"))

## ---- C. harvest and the mortality target ------------------------------------
C1 <- tr[, .(cut_records = sum(STATUSCD == 3, na.rm = TRUE), dead_records = sum(STATUSCD == 2, na.rm = TRUE), alive = sum(STATUSCD == 1, na.rm = TRUE))]
say("[C conus] records on the model plots: alive %d, dead %d, cut (STATUSCD 3) %d", C1$alive, C1$dead, C1$cut_records)
fwrite(C1, file.path(OUT, "C_harvest_conus.csv"))

## ---- D. macroplot --------------------------------------------------------------
D1 <- tr[STATUSCD == 1 & dbh_cm >= PARAMS$min_dbh_cm, .(n = .N, n_macro = sum(geom == "macr"), n_big = sum(dbh_cm >= PARAMS$macro_dbh)),
         by = .(statecd = sub("_.*", "", plotID))][order(-n_macro)]
fwrite(D1, file.path(OUT, "D_macroplot_by_state_conus.csv"))
say("[D conus] macroplot records: %d of %d alive >= 12.7 cm records, on plots of states %s", sum(D1$n_macro), sum(D1$n), paste(D1[n_macro > 0, statecd], collapse = ","))

## ---- E. residuals and rate calibration of the v0.1 holdout -------------------
for (v in names(sets)) {
  s <- sets[[v]]
  M <- fread(file.path("calibration", v, "validation", "predictions.csv"))[set == "holdout"]
  sp <- fread(file.path(DATA_DIR, "prepared", s$prep, s$holdout, "species_dt.csv"))
  ## rate calibration: sums of predicted / observed
  if (v == "oregon") { M[, mort_period_m := mort_m]; M[, plen := s$step] } else M[, plen := period_length]
  E1 <- M[, .(growth_ratio = weighted.mean(growth_m, growth_n, na.rm = TRUE) / weighted.mean(growth, growth_n, na.rm = TRUE),
              deaths_obs = sum(n_died), deaths_pred = sum(mort_period_m * n_at_risk, na.rm = TRUE),
              reg_obs = sum(reg), reg_pred = sum(reg_m),
              ba_obs = sum(ba), ba_pred = sum(ba_m), trees_obs = sum(trees), trees_pred = sum(trees_m), n_at_risk = sum(n_at_risk)),
          by = .(species = species_name)]
  E1[, `:=`(mort_ratio = deaths_pred / deaths_obs, reg_ratio = reg_pred / reg_obs, ba_ratio = ba_pred / ba_obs, trees_ratio = trees_pred / trees_obs)]
  E1 <- E1[order(-n_at_risk)]
  fwrite(E1, file.path(OUT, sprintf("E_rate_calibration_%s.csv", v)))
  tot <- M[, .(growth_ratio = weighted.mean(growth_m, growth_n, na.rm = TRUE) / weighted.mean(growth, growth_n, na.rm = TRUE),
               mort_ratio = sum(mort_period_m * n_at_risk, na.rm = TRUE) / sum(n_died), reg_ratio = sum(reg_m) / sum(reg),
               ba_ratio = sum(ba_m) / sum(ba), trees_ratio = sum(trees_m) / sum(trees))]
  say("[E %s] holdout predicted/observed totals: growth %.2f, deaths %.2f, recruits %.2f, basal area %.3f, stems %.3f",
      v, tot$growth_ratio, tot$mort_ratio, tot$reg_ratio, tot$ba_ratio, tot$trees_ratio)
  ## change metrics against the first inventory
  M[, `:=`(dba_obs = ba - ba_p, dba_pred = ba_m - ba_p, dtr_obs = trees - trees_p, dtr_pred = trees_m - trees_p)]
  E2 <- rbind(cbind(level = "species x plot", variable = "delta basal area", metric(M$dba_obs, M$dba_pred)),
              cbind(level = "species x plot", variable = "delta stems", metric(M$dtr_obs, M$dtr_pred)))
  T <- M[, .(dba_obs = sum(dba_obs), dba_pred = sum(dba_pred), dtr_obs = sum(dtr_obs), dtr_pred = sum(dtr_pred), plen = plen[1]), by = .(siteID, year)]
  E2 <- rbind(E2, cbind(level = "stand total", variable = "delta basal area", metric(T$dba_obs, T$dba_pred)),
              cbind(level = "stand total", variable = "delta stems", metric(T$dtr_obs, T$dtr_pred)))
  E2[, sign_agreement := NA_real_]
  E2[level == "stand total" & variable == "delta basal area", sign_agreement := T[, mean(sign(dba_obs) == sign(dba_pred))]]
  E2[level == "stand total" & variable == "delta stems", sign_agreement := T[, mean(sign(dtr_obs) == sign(dtr_pred))]]
  fwrite(E2, file.path(OUT, sprintf("E_change_metrics_%s.csv", v)))
  say("[E %s] change metrics (holdout): %s", v, paste(sprintf("%s %s R2 %.2f (sign agreement %s)", E2$level, E2$variable, E2$r2,
                                                              ifelse(is.na(E2$sign_agreement), "-", sprintf("%.2f", E2$sign_agreement))), collapse = "; "))
  say("[E %s] stand delta BA per plot: observed mean %.4f, predicted mean %.4f m2/patch; delta stems observed %.3f, predicted %.3f per patch",
      v, mean(T$dba_obs), mean(T$dba_pred), mean(T$dtr_obs), mean(T$dtr_pred))
  ## by interval length (conus) and by initial stand basal area
  if (v == "conus") {
    E3 <- T[, cbind(.(n = .N), metric(dba_obs, dba_pred)[, .(r2, bias)], dba_obs = mean(dba_obs), dba_pred = mean(dba_pred)), by = plen][order(plen)]
    fwrite(E3, file.path(OUT, "E_change_by_interval_conus.csv"))
  }
  S0 <- M[, .(ba0 = sum(ba_p)), by = .(siteID)]
  T <- merge(T, S0, by = "siteID")[, ba0_class := cut(ba0 / PARAMS$patch_size_ha, c(-1, 10, 20, 30, 45, 1e9), labels = c("<10", "10-20", "20-30", "30-45", ">45"))]
  E4 <- T[, cbind(.(n = .N), dba_obs = mean(dba_obs), dba_pred = mean(dba_pred), dtr_obs = mean(dtr_obs), dtr_pred = mean(dtr_pred)), by = ba0_class][order(ba0_class)]
  fwrite(E4, file.path(OUT, sprintf("E_change_by_initial_ba_%s.csv", v)))
  say("[E %s] delta BA by initial stand BA (m2/ha): %s", v, paste(sprintf("%s obs %.4f pred %.4f", E4$ba0_class, E4$dba_obs, E4$dba_pred), collapse = "; "))
  ## per-tree growth by size (conus) and growth by initial dbh class
  if (v == "conus") {
    TR <- fread(file.path("calibration", v, "validation", "tree_predictions.csv"))[set == "holdout" & is.finite(growth)]
    TR[, size := cut(dbh, c(0, 20, 30, 40, 60, 1e9), labels = c("12.7-20", "20-30", "30-40", "40-60", ">60"))]
    E5 <- TR[, cbind(.(n = .N, obs = mean(growth), pred = mean(growth_pred), ratio = mean(growth_pred) / mean(growth)), metric(growth, growth_pred)[, .(r2, spearman)]), by = size][order(size)]
    fwrite(E5, file.path(OUT, "E_tree_growth_by_size_conus.csv"))
    say("[E conus] per-tree growth pred/obs by size: %s", paste(sprintf("%s %.2f (R2 %.2f)", E5$size, E5$ratio, E5$r2), collapse = "; "))
    TRm <- fread(file.path("calibration", v, "validation", "tree_predictions.csv"))[set == "holdout" & is.finite(mort_pred)]
    TRm[, size := cut(dbh, c(0, 20, 30, 40, 60, 1e9), labels = c("12.7-20", "20-30", "30-40", "40-60", ">60"))]
    E6 <- TRm[, .(n = .N, died_obs = mean(died), died_pred = mean(mort_pred), ratio = mean(mort_pred) / mean(died)), by = size][order(size)]
    fwrite(E6, file.path(OUT, "E_tree_mortality_by_size_conus.csv"))
    say("[E conus] per-tree death prob pred/obs by size: %s", paste(sprintf("%s %.2f", E6$size, E6$ratio), collapse = "; "))
  }
  ## recruitment: share of zero cells and the top cells
  E7 <- M[, .(cells = .N, zero_obs = mean(reg == 0), zero_pred = mean(reg_m < 0.5), obs_max = max(reg), pred_max = max(reg_m),
              obs_q99 = quantile(reg, 0.99), pred_q99 = quantile(reg_m, 0.99))]
  fwrite(E7, file.path(OUT, sprintf("E_recruitment_distribution_%s.csv", v)))
  say("[E %s] recruitment cells: observed zero %.2f, predicted < 0.5 %.2f; q99 obs %.1f pred %.1f; max obs %.0f pred %.0f per ha per period",
      v, E7$zero_obs, E7$zero_pred, E7$obs_q99, E7$pred_q99, E7$obs_max, E7$pred_max)
}

## ---- F. convergence, theta, K ------------------------------------------------
for (v in names(sets)) {
  ph <- readRDS(file.path(DATA_DIR, "fits", v, "param_history.rds"))
  lh <- fread(file.path("calibration", v, "loss_history.csv"))
  step <- nrow(lh) / length(ph)
  blocks <- names(ph[[1]])
  F1 <- rbindlist(lapply(blocks, function(b) {
    vals <- lapply(ph, function(s) as.numeric(as.matrix(s[[b]])))
    d <- sapply(seq_along(vals)[-1], function(i) mean(abs(vals[[i]] - vals[[i - 1]])))
    data.table(block = b, n_params = length(vals[[1]]), change_first = d[1], change_mid = d[round(length(d) / 2)], change_last = d[length(d)],
               last_over_first = d[length(d)] / d[1])
  }))
  fwrite(F1, file.path(OUT, sprintf("F_param_change_%s.csv", v)))
  say("[F %s] mean |change| per snapshot (every %g epochs), last / first: %s", v, step,
      paste(sprintf("%s %.2f", F1$block, F1$last_over_first), collapse = "; "))
  lastq <- lh[epoch > 0.8 * max(epoch)]
  say("[F %s] loss over the last fifth of the fit: %s", v, paste(sprintf("%s %.3f -> %.3f", c("dbh", "ba", "trees", "growth", "mort", "reg"),
                                                                       as.numeric(lastq[1, 1:6]), as.numeric(lastq[.N, 1:6])), collapse = "; "))
  p <- torch_load(file.path(DATA_DIR, "fits", v, "parameters.pt"))
  sp <- fread(file.path(DATA_DIR, "prepared", sets[[v]]$prep, sets[[v]]$train, "species_dt.csv"))[order(species)]
  ## FINN's recruitment dispersion: the raw parameter and the theta it implies (softplus, as in finn_class$par_theta_recruits)
  source(file.path("R", "definition.R"))
  m <- build_model(v, p, env_scaling = as.data.frame(fread(file.path("calibration", v, "env_scaling.csv"))))
  theta <- tryCatch(as.numeric(as_array(m$par_theta_recruits)), error = function(e) rep(NA_real_, nrow(sp)))
  K <- exp(as.numeric(as_array(p$reg_logK)))
  F2 <- data.table(species = sp$species_name, theta = theta, K = K)
  if (!is.null(p$reg_adult_logx)) F2[, `:=`(adult_x_cm = exp(as.numeric(as_array(p$reg_adult_logx))), A_half = exp(as.numeric(as_array(p$reg_adult_logA))))]
  fwrite(F2[order(-K)], file.path(OUT, sprintf("F_theta_K_%s.csv", v)))
  say("[F %s] theta: median %.3g, min %.3g, max %.3g; K per ha per step: median %.1f, max %.0f (%s)", v, median(theta, na.rm = TRUE), min(theta, na.rm = TRUE), max(theta, na.rm = TRUE),
      median(K), max(K), F2[which.max(K), species])
}

## ---- G. the "other" class in oregon -----------------------------------------
SB <- fread(file.path("calibration", "oregon", "validation", "bare_ground_succession.csv"))
G1 <- SB[step >= 90, .(ba_ha = mean(ba_ha)), by = .(loc, species_name)][, share := ba_ha / sum(ba_ha), by = loc][order(loc, -share)]
G2 <- G1[, .(top = species_name[1], top_share = share[1], other_share = sum(share[species_name == "other"]), total_ba = sum(ba_ha)), by = loc]
fwrite(G1, file.path(OUT, "G_bare_ground_composition_oregon.csv"))
say("[G oregon] bare-ground end state (last 10 %%): %s", paste(sprintf("%s total %.0f m2/ha, top %s %.0f %%, other %.0f %%", G2$loc, G2$total_ba, G2$top, 100 * G2$top_share, 100 * G2$other_share), collapse = "; "))
Mo <- fread(file.path("calibration", "oregon", "validation", "predictions.csv"))
G3 <- Mo[, .(ba_share_obs = sum(ba) / sum(Mo$ba), reg_share_obs = sum(reg) / sum(Mo$reg), reg_share_pred = sum(reg_m) / sum(Mo$reg_m)), by = species_name][order(-reg_share_pred)][1:6]
fwrite(G3, file.path(OUT, "G_other_shares_oregon.csv"))
say("[G oregon] observed BA share / observed recruit share / predicted recruit share: %s", paste(sprintf("%s %.3f / %.3f / %.3f", G3$species_name, G3$ba_share_obs, G3$reg_share_obs, G3$reg_share_pred), collapse = "; "))

## ---- H. interval spread and tally geometry ----------------------------------
obs <- fread(file.path(DATA_DIR, "prepared", "conus-annual-v1", "train", "obs_dt.csv"))
H1 <- unique(obs[, .(siteID, year, period_length)])[, .N, by = period_length][order(period_length)]
fwrite(H1, file.path(OUT, "H_interval_lengths_conus.csv"))
say("[H conus] intervals (years:n): %s", paste(sprintf("%d:%d", H1$period_length, H1$N), collapse = " "))
writeLines(lines, file.path(OUT, "audit_summary.txt"))
message("done: ", OUT)
