#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/conus/01_fiadb_metadata.R   (variant "conus", step 1 of 9)
##
## WHAT   Reads the entire-US FIADB CSV dump, applies the baseline plot quality
##        filter and writes the per-plot metadata the cohort selection (02) and
##        the FINN preparation (08) work from: one row per plot with its visits,
##        design, composition, demography, condition / management summaries and
##        a management tier; one row per (plot, species) with composition
##        shares; and the species ranking by plot presence.
## WHY    Port of FINN-fia `scripts/01_prepare_metadata.R` (the version that
##        produced the prepared set annual-v5, 2026-09-10), reduced to the
##        tables the chain downstream reads: the diagnostic tables of the
##        original (cohort menu, strictness grid, selection curves) are not
##        written. Every computation that reaches the FINN tables is unchanged.
## INPUT  FIADB CSV dump (PROVENANCE.md): CSV_FIADB_ENTIRE/ENTIRE_{SURVEY,PLOT,
##        SUBPLOT,COND,TREE}.csv, FIADB_REFERENCE/REF_SPECIES.csv, found under
##        <DATA_DIR>/raw/fiadb or the FINN-fia data folder (see FIADB_DIR below).
## OUTPUT <DATA_DIR>/intermediate/conus/plot_metadata.csv          one row per plot
##        <DATA_DIR>/intermediate/conus/plot_observations.csv      one row per (plot, visit)
##        <DATA_DIR>/intermediate/conus/plot_species_summary.csv   one row per (plot, SPCD)
##        <DATA_DIR>/intermediate/conus/species_rank.csv           species by n_plots_present
## RUN    Rscript calibration/conus/01_fiadb_metadata.R   (from the repository root;
##        ~15 min, ~8 GB RAM: the 13 GB TREE file dominates)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))

DATA_DIR <- {
  cands <- c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))
  d <- cands[dir.exists(cands)][1]
  if (is.na(d)) stop("FINN-zoo-data/fia not found (expected two levels above the repository, see calibration/README.md)")
  normalizePath(d)
}
## the FIADB dump: <DATA_DIR>/raw/fiadb/{CSV_FIADB_ENTIRE,FIADB_REFERENCE}, or the
## FINN-fia analysis repository's data folder next to the models folder
FIADB_DIR <- {
  cands <- c(file.path(DATA_DIR, "raw", "fiadb"), file.path("..", "..", "FINN-fia", "data"))
  d <- cands[dir.exists(file.path(cands, "CSV_FIADB_ENTIRE"))][1]
  if (is.na(d)) stop("FIADB dump not found: expected CSV_FIADB_ENTIRE/ under ", paste(cands, collapse = " or "))
  normalizePath(d)
}
CSV_DIR <- file.path(FIADB_DIR, "CSV_FIADB_ENTIRE"); REF_DIR <- file.path(FIADB_DIR, "FIADB_REFERENCE")
OUT <- file.path(DATA_DIR, "intermediate", "conus"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PARAMS <- list(
  ## FIA design geometry (FIA Database User Guide): one subplot = 1/24 acre
  acre_to_ha      = 0.40468564224,
  subplot_area_ha = (1 / 24) * 0.40468564224,
  ## TPA_UNADJ expansion factors that identify the sampling geometry of a tree record
  tpa_subplot     = 6.018046,
  tpa_microplot   = 74.965282,
  tpa_macroplot   = 0.999188,
  ## STATE FIPS codes to names (CONUS + AK + HI)
  statecd = c(1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
              32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56),
  state_abbr = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "DC", "FL", "GA", "HI", "ID", "IL", "IN", "IA", "KS", "KY",
                 "LA", "ME", "MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC", "ND", "OH",
                 "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY"),
  state_name = c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", "Connecticut", "Delaware",
                 "District of Columbia", "Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas",
                 "Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota", "Mississippi",
                 "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire", "New Jersey", "New Mexico", "New York",
                 "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island",
                 "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington",
                 "West Virginia", "Wisconsin", "Wyoming")
)
STATE_LOOKUP <- data.table(STATECD = PARAMS$statecd, STATE_NAME = PARAMS$state_name, STATE_ABBR = PARAMS$state_abbr)
ba_tree_m2 <- function(dia_in) pi * (dia_in * 2.54 / 200)^2

## ---- 1. load ---------------------------------------------------------------
read_fia <- function(file, select = NULL) {
  path <- file.path(CSV_DIR, file)
  if (!file.exists(path)) stop("Missing FIA file: ", path)
  t0 <- Sys.time(); dt <- fread(path, select = select, showProgress = FALSE)
  message(sprintf("[read] %-22s %10d rows (%.0f s)", file, nrow(dt), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  dt
}
survey_dt  <- read_fia("ENTIRE_SURVEY.csv", c("CN", "STATECD", "INVYR", "CYCLE", "SUBCYCLE"))
plot_dt    <- read_fia("ENTIRE_PLOT.csv", c("CN", "SRV_CN", "PREV_PLT_CN", "STATECD", "UNITCD", "COUNTYCD", "PLOT", "INVYR",
                                            "MEASYEAR", "MEASMON", "MEASDAY", "REMPER", "PLOT_STATUS_CD", "QA_STATUS",
                                            "INTENSITY", "DESIGNCD", "LAT", "LON", "ELEV", "MACRO_BREAKPOINT_DIA"))
subplot_dt <- read_fia("ENTIRE_SUBPLOT.csv", c("STATECD", "UNITCD", "COUNTYCD", "PLOT", "SUBP", "INVYR", "SUBP_STATUS_CD",
                                                "MICRCOND", "SUBPCOND", "MACRCOND"))
cond_dt    <- read_fia("ENTIRE_COND.csv", c("STATECD", "UNITCD", "COUNTYCD", "PLOT", "INVYR", "CONDID", "COND_STATUS_CD",
                                            "CONDPROP_UNADJ", "PROP_BASIS", "TRTCD1", "TRTCD2", "TRTCD3", "STDORGCD", "STDAGE",
                                            "DSTRBCD1", "DSTRBCD2", "DSTRBCD3", "RESERVCD", "FORTYPCD", "SLOPE", "ASPECT", "PHYSCLCD"))
ref_species <- fread(file.path(REF_DIR, "REF_SPECIES.csv"), select = c("SPCD", "COMMON_NAME", "GENUS", "SPECIES", "SCIENTIFIC_NAME"),
                     showProgress = FALSE)
tree_dt    <- read_fia("ENTIRE_TREE.csv", c("CN", "PREV_TRE_CN", "STATECD", "UNITCD", "COUNTYCD", "PLOT", "SUBP", "TREE", "INVYR",
                                            "STATUSCD", "AGENTCD", "MORTCD", "STANDING_DEAD_CD", "RECONCILECD", "SPCD", "DIA",
                                            "DIACHECK", "TPA_UNADJ"))

## ---- 2. baseline filter and keys -------------------------------------------
mkPLOT <- function(dt) dt[, plotID := paste(STATECD, UNITCD, COUNTYCD, PLOT, sep = "_")]
mkSUBP <- function(dt) dt[, subplotID := paste(STATECD, UNITCD, COUNTYCD, PLOT, SUBP, sep = "_")]
mkPLOT(plot_dt); mkPLOT(subplot_dt); mkSUBP(subplot_dt); mkPLOT(cond_dt); mkPLOT(tree_dt); mkSUBP(tree_dt)

valid_srv <- survey_dt[INVYR != 9999, CN]
plot_sel  <- plot_dt[SRV_CN %in% valid_srv & PLOT_STATUS_CD == 1 & QA_STATUS == 1 & INTENSITY == 1 & INVYR != 9999 &
                     !is.na(LAT) & !is.na(LON)]
message(sprintf("[filter] %d / %d PLOT rows kept (%d unique plots)", nrow(plot_sel), nrow(plot_dt), uniqueN(plot_sel$plotID)))
keep <- unique(plot_sel[, .(plotID, INVYR)]); setkey(keep, plotID, INVYR)
subplot_sel <- subplot_dt[keep, on = .(plotID, INVYR), nomatch = NULL]
cond_sel    <- cond_dt[keep, on = .(plotID, INVYR), nomatch = NULL]
tree_sel    <- tree_dt[keep, on = .(plotID, INVYR), nomatch = NULL]
rm(plot_dt, subplot_dt, cond_dt, tree_dt); invisible(gc(verbose = FALSE))
tree_sel[, geom := fifelse(abs(TPA_UNADJ - PARAMS$tpa_subplot) < 0.01, "subp",
                  fifelse(abs(TPA_UNADJ - PARAMS$tpa_microplot) < 0.5, "micr",
                  fifelse(abs(TPA_UNADJ - PARAMS$tpa_macroplot) < 0.01, "macr", "other")))]

## ---- 3. per-(plot, visit) tree-side aggregates -----------------------------
subp_per_v <- subplot_sel[SUBP_STATUS_CD == 1, .(n_subp_sampled = uniqueN(SUBP),
                                                 subplots_present = paste(sort(unique(SUBP)), collapse = ",")), by = .(plotID, INVYR)]
tree_sel[, ba_m2 := ba_tree_m2(DIA)]
live_pvs <- tree_sel[geom == "subp" & STATUSCD == 1 & !is.na(SPCD), .(n_alive = .N, ba_m2 = sum(ba_m2, na.rm = TRUE)), by = .(plotID, INVYR, SPCD)]
live_pv  <- live_pvs[, .(n_alive_total = sum(n_alive), ba_m2_total = sum(ba_m2)), by = .(plotID, INVYR)]
live_pv  <- merge(live_pv, subp_per_v, by = c("plotID", "INVYR"), all.x = TRUE)
live_pv[is.na(n_subp_sampled), n_subp_sampled := 0L]
live_pv[, area_ha_subp := n_subp_sampled * PARAMS$subplot_area_ha]
live_pv[, ba_m2_ha := fifelse(area_ha_subp > 0, ba_m2_total / area_ha_subp, NA_real_)]
live_pvs <- merge(live_pvs, live_pv[, .(plotID, INVYR, ba_m2_total, n_alive_total)], by = c("plotID", "INVYR"))
live_pvs[, ba_share := fifelse(ba_m2_total > 0, ba_m2 / ba_m2_total, NA_real_)]
live_pvs[, count_share := fifelse(n_alive_total > 0, n_alive / n_alive_total, NA_real_)]
setorder(live_pvs, plotID, INVYR, -ba_share); live_pvs[, rank_ba := seq_len(.N), by = .(plotID, INVYR)]
setorder(live_pvs, plotID, INVYR, -count_share); live_pvs[, rank_count := seq_len(.N), by = .(plotID, INVYR)]
setorder(live_pvs, plotID, INVYR, rank_ba)

tree_chron <- tree_sel[geom == "subp" & !is.na(SPCD)]
setorder(tree_chron, plotID, subplotID, TREE, INVYR)
tree_chron[, `:=`(statuscd_prev = shift(STATUSCD, 1L, type = "lag"), invyr_prev = shift(INVYR, 1L, type = "lag")), by = .(subplotID, TREE)]
tree_chron[, died := !is.na(statuscd_prev) & statuscd_prev == 1 & STATUSCD == 2]
tree_chron[, died_silvi := died & !is.na(AGENTCD) & AGENTCD == 80]
tree_chron[, died_natural := died & !died_silvi]
tree_chron[, cut_event := STATUSCD == 3]
tree_chron[, record_idx := seq_len(.N), by = .(subplotID, TREE)]
tree_chron[, first_invyr_plot := min(INVYR), by = plotID]
tree_chron[, recruited := record_idx == 1 & INVYR > first_invyr_plot & STATUSCD == 1]
tree_chron[!is.na(RECONCILECD) & RECONCILECD %in% c(6, 7, 8), recruited := FALSE]
demog_pv <- tree_chron[, .(n_mort_total = sum(died, na.rm = TRUE), n_mort_natural = sum(died_natural, na.rm = TRUE),
                           n_mort_silvi = sum(died_silvi, na.rm = TRUE), n_cut_events = sum(cut_event, na.rm = TRUE),
                           n_recruit = sum(recruited, na.rm = TRUE),
                           n_diacheck_ok = sum(STATUSCD == 1 & !is.na(DIACHECK) & DIACHECK == 0),
                           n_live_with_dc = sum(STATUSCD == 1 & !is.na(DIACHECK))), by = .(plotID, INVYR)]
rm(tree_chron); invisible(gc(verbose = FALSE))

## ---- 4. per-(plot, visit) condition-side aggregates ------------------------
for (col in c("TRTCD1", "TRTCD2", "TRTCD3", "DSTRBCD1", "DSTRBCD2", "DSTRBCD3")) cond_sel[is.na(get(col)), (col) := 0]
setorder(cond_sel, plotID, INVYR, -CONDPROP_UNADJ)
dom_cond <- cond_sel[, .SD[1L], by = .(plotID, INVYR),
                     .SDcols = c("CONDID", "COND_STATUS_CD", "CONDPROP_UNADJ", "PROP_BASIS", "STDORGCD", "STDAGE", "FORTYPCD",
                                 "RESERVCD", "SLOPE", "ASPECT", "PHYSCLCD")]
setnames(dom_cond, c("CONDID", "COND_STATUS_CD", "CONDPROP_UNADJ", "PROP_BASIS", "STDORGCD", "STDAGE", "FORTYPCD", "RESERVCD",
                     "SLOPE", "ASPECT", "PHYSCLCD"),
         c("dom_condid", "dom_cond_status", "dom_condprop", "prop_basis_visit", "stdorgcd", "stdage", "fortypcd", "reservcd",
           "slope", "aspect", "physclcd"))
n_cond_pv <- cond_sel[, .(n_conditions = uniqueN(CONDID)), by = .(plotID, INVYR)]
mgmt_pv <- cond_sel[, .(trtcd_any = any(TRTCD1 != 0 | TRTCD2 != 0 | TRTCD3 != 0),
                        trtcd_cut = any(TRTCD1 == 10 | TRTCD2 == 10 | TRTCD3 == 10),
                        dstrb_any = any(DSTRBCD1 != 0 | DSTRBCD2 != 0 | DSTRBCD3 != 0),
                        dstrb_codes = paste(sort(unique(c(DSTRBCD1, DSTRBCD2, DSTRBCD3))), collapse = ",")), by = .(plotID, INVYR)]
mgmt_pv[, dstrb_codes := gsub("(^|,)0(,|$)", "\\2", dstrb_codes)]
mgmt_pv[dstrb_codes == "", dstrb_codes := NA_character_]
cond_pv <- Reduce(function(a, b) merge(a, b, by = c("plotID", "INVYR"), all = TRUE), list(dom_cond, n_cond_pv, mgmt_pv))

## ---- 5. plot_observations.csv ----------------------------------------------
plot_visit <- plot_sel[, .(plotID, STATECD, INVYR, MEASYEAR, REMPER, LAT, LON, ELEV, DESIGNCD, PREV_PLT_CN)]
top1_pv <- live_pvs[rank_ba == 1L, .(plotID, INVYR, top1_spp_at_visit = SPCD, top1_share_at_visit = round(ba_share, 4))]
plot_obs <- Reduce(function(a, b) merge(a, b, by = c("plotID", "INVYR"), all.x = TRUE), list(plot_visit, live_pv, demog_pv, top1_pv, cond_pv))
setorder(plot_obs, plotID, MEASYEAR)
plot_obs[, observation_id := seq_len(.N), by = plotID]
plot_obs[, remper_since_prev := REMPER]
plot_obs[observation_id == 1L, remper_since_prev := NA_real_]
obs_cols <- intersect(c("plotID", "observation_id", "INVYR", "MEASYEAR", "remper_since_prev", "n_subp_sampled", "area_ha_subp",
                        "n_alive_total", "ba_m2_total", "ba_m2_ha", "n_mort_total", "n_mort_natural", "n_mort_silvi", "n_cut_events",
                        "n_recruit", "n_diacheck_ok", "n_live_with_dc", "top1_spp_at_visit", "top1_share_at_visit", "dom_condid",
                        "dom_cond_status", "dom_condprop", "prop_basis_visit", "n_conditions", "stdorgcd", "stdage", "fortypcd",
                        "reservcd", "slope", "aspect", "physclcd", "trtcd_any", "trtcd_cut", "dstrb_any", "dstrb_codes", "STATECD",
                        "DESIGNCD", "LAT", "LON", "ELEV", "PREV_PLT_CN"), names(plot_obs))
plot_obs <- plot_obs[, ..obs_cols]
fwrite(plot_obs, file.path(OUT, "plot_observations.csv"))

## ---- 6. plot_species_summary.csv -------------------------------------------
n_visits <- plot_obs[, .(n_visits = .N), by = plotID]
sp_sum <- live_pvs[, .(ba_share_sum = sum(ba_share, na.rm = TRUE), count_share_sum = sum(count_share, na.rm = TRUE),
                       ba_m2_sum = sum(ba_m2, na.rm = TRUE), n_alive_sum = sum(n_alive, na.rm = TRUE), n_visits_alive = .N), by = .(plotID, SPCD)]
sp_avg <- merge(sp_sum, n_visits, by = "plotID")
sp_avg[, `:=`(ba_share_avg = ba_share_sum / n_visits, count_share_avg = count_share_sum / n_visits,
              ba_m2_avg = ba_m2_sum / n_visits, n_alive_avg = n_alive_sum / n_visits)]
sp_avg[, c("ba_share_sum", "count_share_sum", "ba_m2_sum", "n_alive_sum", "n_visits") := NULL]
first_visit_invyr <- plot_obs[, .SD[which.min(MEASYEAR)], by = plotID, .SDcols = "INVYR"]
sp_first <- live_pvs[first_visit_invyr, on = c("plotID", "INVYR"), nomatch = NULL][
  , .(plotID, SPCD, ba_share_first = ba_share, count_share_first = count_share, ba_m2_first = ba_m2, n_trees_first = n_alive, rank_ba_first = rank_ba)]
sp_long <- merge(sp_avg, sp_first, by = c("plotID", "SPCD"), all.x = TRUE)
sp_long[is.na(ba_share_first), ba_share_first := 0]; sp_long[is.na(count_share_first), count_share_first := 0]
sp_long[is.na(ba_m2_first), ba_m2_first := 0]; sp_long[is.na(n_trees_first), n_trees_first := 0L]
setorder(sp_long, plotID, -ba_share_avg); sp_long[, rank_ba_avg := seq_len(.N), by = plotID]
setorder(sp_long, plotID, -count_share_avg); sp_long[, rank_count_avg := seq_len(.N), by = plotID]
setorder(sp_long, plotID, rank_ba_avg)
plot_species_summary <- merge(
  sp_long[, .(plotID, SPCD, ever_alive = 1L, n_visits_alive, ba_share_first = round(ba_share_first, 4), ba_share_avg = round(ba_share_avg, 4),
              count_share_first = round(count_share_first, 4), count_share_avg = round(count_share_avg, 4),
              ba_m2_first = round(ba_m2_first, 3), ba_m2_avg = round(ba_m2_avg, 3), n_trees_first, n_trees_avg = round(n_alive_avg, 2),
              rank_ba_first, rank_ba_avg)],
  ref_species[, .(SPCD, SCIENTIFIC_NAME, COMMON_NAME)], by = "SPCD", all.x = TRUE)
setcolorder(plot_species_summary, c("plotID", "SPCD", "SCIENTIFIC_NAME", "COMMON_NAME"))
setorder(plot_species_summary, plotID, rank_ba_avg)
fwrite(plot_species_summary, file.path(OUT, "plot_species_summary.csv"))

## ---- 7. plot_metadata.csv (one row per plot) -------------------------------
id_pp <- plot_obs[, .SD[which.min(MEASYEAR)], by = plotID, .SDcols = c("STATECD", "LAT", "LON", "ELEV")]
id_pp <- merge(id_pp, STATE_LOOKUP, by = "STATECD", all.x = TRUE)
setnames(id_pp, c("LAT", "LON", "ELEV"), c("lat", "lon", "elev_m"))
tmp_pp <- plot_obs[, {
  yrs <- sort(unique(MEASYEAR)); diffs <- if (length(yrs) > 1) diff(yrs) else NA_real_
  list(n_inventories_raw = length(yrs), first_year = min(yrs), last_year = max(yrs), time_span_years = max(yrs) - min(yrs),
       mean_remper = if (length(yrs) > 1) mean(diffs, na.rm = TRUE) else NA_real_,
       sd_remper = if (length(yrs) > 2) sd(diffs, na.rm = TRUE) else NA_real_,
       inv_year_list = paste(yrs, collapse = ","), record_begins_pre2000 = as.integer(min(yrs) < 2000))
}, by = plotID]
tmp_pp[, panel_cycle := fifelse(n_inventories_raw == 1L, "single",
                       fifelse(mean_remper >= 4 & mean_remper <= 6, "5yr",
                       fifelse(mean_remper >= 8 & mean_remper <= 12, "10yr", "irregular")))]
tmp_pp[, is_consecutive_panel := as.integer(is.na(sd_remper) | sd_remper <= 1.5)]
design_pp <- tree_sel[, .(has_subp_dbh5p_0p067ha = as.integer(any(geom == "subp")), has_micr_dbh1to5_0p005ha = as.integer(any(geom == "micr")),
                          has_macr_dbh20p_0p405ha = as.integer(any(geom == "macr"))), by = plotID]
nsubp_first <- plot_obs[order(MEASYEAR), .SD[1L], by = plotID, .SDcols = c("n_subp_sampled", "area_ha_subp")]
setnames(nsubp_first, c("n_subp_sampled", "area_ha_subp"), c("n_subp_sampled_first", "area_ha_subp_first"))
design_pp <- merge(design_pp, nsubp_first, by = "plotID", all.x = TRUE)
designcd_pp <- plot_obs[, .(designcd_first = DESIGNCD[which.min(MEASYEAR)], designcd_constant_across_visits = as.integer(uniqueN(DESIGNCD) == 1L)), by = plotID]
design_pp <- merge(design_pp, designcd_pp, by = "plotID", all.x = TRUE)
spc_pp <- subp_per_v[, .(subplots_consistent_over_time = as.integer(uniqueN(subplots_present) == 1L)), by = plotID]
design_pp <- merge(design_pp, spc_pp, by = "plotID", all.x = TRUE)

tba_first <- sp_long[rank_ba_first <= 3L & !is.na(rank_ba_first), .(plotID, rank_ba_first, SPCD, ba_share_first)]
tba_first_w <- dcast(tba_first, plotID ~ rank_ba_first, value.var = c("SPCD", "ba_share_first"))
setnames(tba_first_w, old = c("SPCD_1", "SPCD_2", "SPCD_3", "ba_share_first_1", "ba_share_first_2", "ba_share_first_3"),
         new = c("top1_species_ba", "top2_species_ba", "top3_species_ba", "top1_share_ba_first", "top2_share_ba_first", "top3_share_ba_first"),
         skip_absent = TRUE)
tba_avg <- sp_long[rank_ba_avg <= 3L, .(plotID, rank_ba_avg, SPCD, ba_share_avg)]
tba_avg_w <- dcast(tba_avg, plotID ~ rank_ba_avg, value.var = c("SPCD", "ba_share_avg"))
setnames(tba_avg_w, old = c("ba_share_avg_1", "ba_share_avg_2", "ba_share_avg_3"), new = c("top1_share_ba_avg", "top2_share_ba_avg", "top3_share_ba_avg"), skip_absent = TRUE)
tba_avg_w[, c("SPCD_1", "SPCD_2", "SPCD_3") := NULL]
tcn_avg <- sp_long[rank_count_avg <= 3L, .(plotID, rank_count_avg, SPCD, count_share_avg)]
tcn_avg_w <- dcast(tcn_avg, plotID ~ rank_count_avg, value.var = c("SPCD", "count_share_avg"))
setnames(tcn_avg_w, old = c("SPCD_1", "SPCD_2", "SPCD_3", "count_share_avg_1", "count_share_avg_2", "count_share_avg_3"),
         new = c("top1_species_count", "top2_species_count", "top3_species_count", "top1_share_count_avg", "top2_share_count_avg", "top3_share_count_avg"),
         skip_absent = TRUE)
sp_basics <- sp_long[, .(n_species_ever = .N, species_list_alive = paste(sort(unique(SPCD)), collapse = ",")), by = plotID]
ba_totals <- plot_obs[, .(ba_total_first_m2_ha = ba_m2_ha[which.min(MEASYEAR)], ba_total_avg_m2_ha = mean(ba_m2_ha, na.rm = TRUE)), by = plotID]
comp_pp <- Reduce(function(a, b) merge(a, b, by = "plotID", all = TRUE), list(tba_first_w, tba_avg_w, tcn_avg_w, sp_basics, ba_totals))

demog_pp <- plot_obs[, .(n_alive_first = n_alive_total[which.min(MEASYEAR)], n_alive_last = n_alive_total[which.max(MEASYEAR)],
                         n_mort_total_window = sum(n_mort_total, na.rm = TRUE), n_mort_natural_window = sum(n_mort_natural, na.rm = TRUE),
                         n_mort_silvi_window = sum(n_mort_silvi, na.rm = TRUE), n_recruit_window = sum(n_recruit, na.rm = TRUE)), by = plotID]
demog_pp <- merge(demog_pp, tmp_pp[, .(plotID, time_span_years)], by = "plotID")
demog_pp <- merge(demog_pp, design_pp[, .(plotID, area_ha_subp_first)], by = "plotID", all.x = TRUE)
demog_pp[, mortality_rate_natural_pct_yr := fifelse(n_alive_first > 0 & time_span_years > 0, 100 * n_mort_natural_window / n_alive_first / time_span_years, NA_real_)]
demog_pp[, recruitment_rate_ha_yr := fifelse(!is.na(area_ha_subp_first) & area_ha_subp_first > 0 & time_span_years > 0,
                                             n_recruit_window / area_ha_subp_first / time_span_years, NA_real_)]
demog_pp[, c("time_span_years", "area_ha_subp_first") := NULL]

parse_dstrb_bucket <- function(s, range) {
  if (is.na(s) || s == "") return(0L)
  as.integer(any(as.integer(strsplit(s, ",")[[1]]) %in% range, na.rm = TRUE))
}
site_pp <- plot_obs[order(MEASYEAR), .(
  forest_type_first = fortypcd[1L], forest_type_last = fortypcd[.N], stdage_first = stdage[1L], stdage_last = stdage[.N],
  stdorgcd_first = stdorgcd[1L], stdorgcd_any_artificial = as.integer(any(stdorgcd == 1L, na.rm = TRUE)),
  slope_pct = mean(slope, na.rm = TRUE), aspect_deg = mean(aspect, na.rm = TRUE), physclcd = physclcd[1L],
  reservcd_first = reservcd[1L], is_reserved_ever = as.integer(any(reservcd == 1L, na.rm = TRUE)),
  trtcd_any_recorded = as.integer(any(trtcd_any)), trtcd_cut_recorded = as.integer(any(trtcd_cut)),
  n_conditions_first = n_conditions[1L], dominant_condprop_first = dom_condprop[1L],
  single_condition_plot = as.integer(n_conditions[1L] == 1L), cond_homogeneous_95 = as.integer(dom_condprop[1L] >= 0.95),
  had_disturbance = as.integer(any(dstrb_any)),
  disturbance_codes = { codes <- unique(unlist(strsplit(stats::na.omit(dstrb_codes), ","))); codes <- codes[codes != "" & codes != "0"]
                        if (length(codes) > 0) paste(sort(codes), collapse = ",") else NA_character_ },
  prop_basis = paste(sort(unique(prop_basis_visit)), collapse = ",")), by = plotID]
site_pp[, stand_age_tier_first := fcase(is.na(stdage_first), NA_character_, stdage_first < 30, "young", stdage_first < 80, "mid",
                                        stdage_first < 150, "mature", default = "old_growth")]
site_pp[, had_fire := vapply(disturbance_codes, function(x) parse_dstrb_bucket(x, c(30, 31, 32)), integer(1))]
site_pp[, had_insect_disease := vapply(disturbance_codes, function(x) parse_dstrb_bucket(x, c(10:22)), integer(1))]
site_pp[, had_weather := vapply(disturbance_codes, function(x) parse_dstrb_bucket(x, c(50:54)), integer(1))]
qf_pp <- plot_obs[, .(pct_trees_diacheck_ok = sum(n_diacheck_ok, na.rm = TRUE) / pmax(sum(n_live_with_dc, na.rm = TRUE), 1),
                      pct_subplots_fully_sampled = mean(n_subp_sampled, na.rm = TRUE) / 4), by = plotID]
prev_chain_pp <- plot_sel[order(plotID, MEASYEAR), .(prev_plt_cn_chain_valid = as.integer(all(!is.na(PREV_PLT_CN[-1L])) | .N == 1L)), by = plotID]
first_v <- plot_obs[order(MEASYEAR), .SD[1L], by = plotID, .SDcols = c("dom_cond_status", "n_alive_total")]
setnames(first_v, c("dom_cond_status", "n_alive_total"), c("dom_cs_first", "n_alive_first_check"))
empty_pp <- merge(first_v, site_pp[, .(plotID, had_disturbance, trtcd_cut_recorded)], by = "plotID", all.x = TRUE)
empty_pp[, is_empty_forest_first := as.integer(!is.na(dom_cs_first) & dom_cs_first == 1L & (is.na(n_alive_first_check) | n_alive_first_check == 0L))]
empty_pp[, likely_post_disturbance := as.integer(is_empty_forest_first == 1L & had_disturbance == 1L)]
empty_pp[, likely_post_clearcut := as.integer(is_empty_forest_first == 1L & trtcd_cut_recorded == 1L)]
empty_pp <- empty_pp[, .(plotID, is_empty_forest_first, likely_post_disturbance, likely_post_clearcut)]
cut_pp   <- tree_sel[STATUSCD == 3, .(n_cut_trees_total = .N), by = plotID]
silvi_pp <- tree_sel[STATUSCD == 2 & !is.na(AGENTCD) & AGENTCD == 80, .(n_silvi_mort_total = .N), by = plotID]
mgmt_tier <- merge(site_pp[, .(plotID, stdorgcd_first, trtcd_cut_recorded, trtcd_any_recorded)], demog_pp[, .(plotID, n_mort_silvi_window)], by = "plotID")
mgmt_tier <- merge(mgmt_tier, cut_pp, by = "plotID", all.x = TRUE)
mgmt_tier <- merge(mgmt_tier, silvi_pp, by = "plotID", all.x = TRUE)
mgmt_tier[is.na(n_cut_trees_total), n_cut_trees_total := 0L]; mgmt_tier[is.na(n_silvi_mort_total), n_silvi_mort_total := 0L]
mgmt_tier[, management_tier := fcase(
  stdorgcd_first == 0 & trtcd_any_recorded == 0 & n_cut_trees_total == 0 & n_silvi_mort_total == 0, "M0_pristine",
  stdorgcd_first == 0 & trtcd_cut_recorded == 0 & n_cut_trees_total == 0 & n_silvi_mort_total == 0, "M1_natural_origin",
  trtcd_cut_recorded == 0 & n_cut_trees_total == 0 & n_silvi_mort_total == 0, "M2_no_cut",
  default = "M3_managed")]
md <- Reduce(function(a, b) merge(a, b, by = "plotID", all.x = TRUE),
             list(id_pp, tmp_pp, design_pp, comp_pp, demog_pp, site_pp, qf_pp, prev_chain_pp, empty_pp, mgmt_tier))
md[, n_inventories_clean := n_inventories_raw]
md[, removal_reason := NA_character_]
fwrite(md, file.path(OUT, "plot_metadata.csv"))
message(sprintf("[out] plot_metadata.csv: %d plots x %d columns", nrow(md), ncol(md)))

## ---- 8. species ranking by plot presence -----------------------------------
species_rank <- sp_long[, .(n_plots_present = uniqueN(plotID), total_ba_share_sum = round(sum(ba_share_avg, na.rm = TRUE), 4),
                            mean_ba_share_when_present = round(mean(ba_share_avg, na.rm = TRUE), 4)), by = SPCD]
setorder(species_rank, -n_plots_present)
species_rank[, rank_global := seq_len(.N)]
species_rank <- merge(species_rank, ref_species[, .(SPCD, SCIENTIFIC_NAME, COMMON_NAME)], by = "SPCD", all.x = TRUE)
setcolorder(species_rank, c("rank_global", "SPCD", "SCIENTIFIC_NAME", "COMMON_NAME", "n_plots_present", "total_ba_share_sum", "mean_ba_share_when_present"))
setorder(species_rank, rank_global)
fwrite(species_rank, file.path(OUT, "species_rank.csv"))
message(sprintf("[out] species_rank.csv: %d species; plot_species_summary.csv: %d rows; plot_observations.csv: %d rows",
                nrow(species_rank), nrow(plot_species_summary), nrow(plot_obs)))
