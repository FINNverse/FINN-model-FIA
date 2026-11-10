#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/prepare_conus_v2_annual.R   (variant "conus", data v2, step A)
##
## V0.2 DATA (2026-11-10): copy of calibration/conus/08_prepare_annual_v5.R with
##        ONE change, the recruit rule (audit finding B). A tree that was tallied
##        on the microplot below 12.7 cm at the previous visit carries a
##        PREV_TRE_CN; the sapling filter drops that earlier record, and v0.1
##        then called the tree "alive" at the first census it appears >= 12.7 cm,
##        so it never counted as ingrowth (720 of 8723 ingrowth trees on the
##        model plots, 8 %). Here a living tree is "new" at the FIRST RETAINED
##        record of its chronology (whether or not a sapling record preceded
##        it), "alive" afterwards; so `reg` is the complete ingrowth into the
##        12.7 cm tally and the closed mortality cohort is unchanged. Output
##        version annual-v6-ingrowth. Everything else is the step-8 header below.
##
## WHAT   Turns the cohort tree records (step 3), the static WorldClim table
##        (step 4), the yearly TerraClimate table (step 6) and the species
##        ranking (step 7) into the FINN input tables of the prepared set
##        `annual-v5`: the 20 top-ranked species, plots with four consistently
##        sampled subplots whose 20 species hold >= 99 % of the first-inventory
##        basal area, macroplot thinning, trees >= 12.7 cm, FINN::makeObsData()
##        per subplot and FINN::resolveSiteIDs(), an annual time axis with each
##        plot's real remeasurement intervals (period_length), annualised growth
##        and mortality, a tree-level observation table, yearly climate per
##        site, 2-degree spatial blocks with 20 % held out, and the three sets
##        national / train1500 / holdout_blocked.
## WHY    Port of FINN-fia `scripts/07_prepare_finn_inputs.R` at the commit
##        that wrote annual-v5 (58fcf3e, 2026-09-10), same PARAMS, same
##        computations, same seeds, so the ported chain can be checked against
##        the FINN-fia tables (calibration/expected_preparation.csv). Only the
##        paths changed, and the size-class / quantile columns that annual-v5
##        carried for other experiments are still written (they are ignored by
##        the fit). The sapling-free, land-use-free, complete-site logic is the
##        original's.
## INPUT  <DATA_DIR>/intermediate/conus/{calibration_set_plots.csv, calibration_set_species_v.csv,
##        calibration_set_trees.csv, plot_env.csv, annual_climate_by_plot.csv, selected_species.csv}
## OUTPUT <DATA_DIR>/prepared/annual-v5/{national,train1500,holdout_blocked}/
##          obs_dt.csv, obs_dt_patches.csv, tree_obs_dt.csv, env_dt.csv (z-scaled), env_unscaled_dt.csv,
##          env_scales.csv, init_cohort_dt.csv, species_dt.csv, siteID_dt.csv, site_calendar.csv
##        <DATA_DIR>/prepared/annual-v5/{PREPARATION.md, site_blocks.csv, holdout_separation.csv}
## RUN    Rscript calibration/conus/08_prepare_annual_v5.R   (from the repository root, ~5 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(FINN) })

DATA_DIR <- {
  cands <- c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))
  d <- cands[dir.exists(cands)][1]
  if (is.na(d)) stop("FINN-zoo-data/fia not found (expected two levels above the repository, see calibration/README.md)")
  normalizePath(d)
}
INT <- file.path(DATA_DIR, "intermediate", "conus")

PARAMS <- list(
  max_species           = 20L,      # the top 20 of the step-7 ranking: strict, so the 99 % rule admits only plots made of modelled species
  min_selected_ba_share = 0.99,     # a plot is kept if the 20 species hold >= 99 % of its first-inventory basal area
  n_patches             = 4L,       # the four subplots of a plot
  patch_size_ha         = 0.404686 / 6.018046 / 4,   # one subplot = a quarter of the 1/6.018-acre cluster = 0.01681 ha (FINN-fia's acre constant, so `reg` reproduces to the last digit)
  min_n_years           = 2L,       # makeObsData: >= 2 inventories with equal coverage on all patches
  fix_period_length     = NULL,     # mixed remeasurement intervals kept (per-plot period_length)
  dbh_growth_thresh     = c(-10, 50),   # drop plots with an impossible dbh change (cm)
  env_vars              = c("tavg_c", "prec_mm"),   # static WorldClim normals: used to drop plots outside a cell
  xy_vars               = c("x", "y"),              # plot coordinates, carried in env for the blocks and maps
  min_dbh_cm            = 12.7,     # the subplot tally threshold; saplings dropped
  annual_vars           = c("tavg_c", "prec_mm", "cwd_z"),   # yearly TerraClimate climate per plot (step 6)
  macroplot_keep_one_in = 6L,       # Pacific-Northwest macroplot records represent 1/6 of the subplot density
  max_span_years        = 25L,      # drop plots whose first-to-last inventory span exceeds the horizon
  size_class_breaks     = c(20, 30, 40, 60),   # diagnostic stem-count classes carried on obs_dt (not fitted)
  dbh_probs             = c(0.1, 0.5, 0.9),    # diagnostic diameter quantiles carried on obs_dt (not fitted)
  prep_version          = "annual-v6-ingrowth",
  block_deg             = 2,        # spatial blocks of 2 x 2 degrees
  holdout_block_share   = 0.2,      # 20 % of the blocks held out
  block_seed            = 21L,
  sets                  = list(national        = list(),
                               train1500       = list(n = 1500L, seed = 1L,  blocks = "train"),
                               holdout_blocked = list(n = 1000L, seed = 11L, blocks = "holdout"))
)
OUT_PREP <- file.path(DATA_DIR, "prepared", PARAMS$prep_version)
dir.create(OUT_PREP, recursive = TRUE, showWarnings = FALSE)
ENV_COLS <- c(PARAMS$env_vars, PARAMS$xy_vars)

thin_macroplot_trees <- function(dt, keep_one_in) {
  is_macr <- dt$geom == "macr"
  if (!any(is_macr, na.rm = TRUE)) return(dt)
  key <- paste(dt$plotID, dt$SUBP, dt$TREE, sep = "_")[which(is_macr)]
  h <- vapply(unique(key), function(s) { v <- utf8ToInt(s); sum(v * seq_along(v)) %% keep_one_in }, numeric(1))
  keep <- h[key] == 0
  drop_idx <- which(is_macr)[!keep]
  message(sprintf("[macroplot] %d macroplot records on %d plots; keeping 1 in %d per tree -> dropped %d, kept %d",
                  sum(is_macr, na.rm = TRUE), uniqueN(dt$plotID[which(is_macr)]), keep_one_in, length(drop_idx),
                  sum(is_macr, na.rm = TRUE) - length(drop_idx)))
  if (length(drop_idx)) dt <- dt[-drop_idx]
  dt
}

## ---- 1. inputs -------------------------------------------------------------
plots_dt   <- fread(file.path(INT, "calibration_set_plots.csv"))
species_dt <- fread(file.path(INT, "calibration_set_species_v.csv"))
env_in_dt  <- fread(file.path(INT, "plot_env.csv"))
message(sprintf("[load] %d plots, %d species in V, env for %d plots", nrow(plots_dt), nrow(species_dt), nrow(env_in_dt)))

## ---- 2. four consistently sampled subplots ---------------------------------
n_total <- nrow(plots_dt)
plots_dt <- plots_dt[n_subp_sampled_first == 4L & subplots_consistent_over_time == TRUE]
message(sprintf("[filter] all 4 subplots usable: %d -> %d plots", n_total, nrow(plots_dt)))
pilot_plotIDs <- plots_dt$plotID

## ---- 3. the species: top 20 of the step-7 ranking --------------------------
sel_file <- fread(file.path(INT, "selected_species.csv"))[order(rank)][seq_len(min(PARAMS$max_species, .N))]
selected_spcds <- unique(sel_file$SPCD)
selected_species <- unique(species_dt[SPCD %in% selected_spcds, .(SPCD, SCIENTIFIC_NAME)])
message(sprintf("[species] %d species: %s", length(selected_spcds), paste(sel_file$SCIENTIFIC_NAME, collapse = ", ")))

## ---- 4. trees of the plots ---------------------------------------------------
tree_cols <- c("plotID", "subplotID", "STATECD", "COUNTYCD", "PLOT", "SUBP", "TREE", "CN", "PREV_TRE_CN", "INVYR", "SPCD",
               "SCIENTIFIC_NAME", "STATUSCD", "AGENTCD", "MORTCD", "RECONCILECD", "dbh_cm", "ba_m2", "height_m", "TPA_UNADJ", "geom")
trees_raw <- fread(file.path(INT, "calibration_set_trees.csv"), select = tree_cols)[plotID %in% pilot_plotIDs]
trees_raw <- thin_macroplot_trees(trees_raw, PARAMS$macroplot_keep_one_in)
n_before <- nrow(trees_raw)
trees_raw <- trees_raw[is.na(dbh_cm) | dbh_cm >= PARAMS$min_dbh_cm]
message(sprintf("[filter] dbh >= %.1f cm: %d -> %d records on %d plots", PARAMS$min_dbh_cm, n_before, nrow(trees_raw), uniqueN(trees_raw$plotID)))

## ---- 5. the 99 % purity rule at the first inventory ------------------------
first_yr <- trees_raw[, .(first_yr = min(INVYR)), by = plotID]
ba_first <- merge(trees_raw[STATUSCD == 1L, .(plotID, INVYR, SPCD, ba_m2)], first_yr, by = "plotID")[INVYR == first_yr]
site_share <- ba_first[, .(ba_total = sum(ba_m2, na.rm = TRUE), ba_selected = sum(ba_m2[SPCD %in% selected_spcds], na.rm = TRUE)), by = plotID][
  , share_selected := fifelse(ba_total > 0, ba_selected / ba_total, 0)]
keep_plotIDs <- site_share[share_selected >= PARAMS$min_selected_ba_share, plotID]
message(sprintf("[filter] %.0f %% purity: %d -> %d plots", 100 * PARAMS$min_selected_ba_share, length(pilot_plotIDs), length(keep_plotIDs)))
trees_raw <- trees_raw[plotID %in% keep_plotIDs]
n_before <- nrow(trees_raw)
trees_raw <- trees_raw[SPCD %in% selected_spcds]
message(sprintf("[filter] residual non-selected trees dropped: %d -> %d records", n_before, nrow(trees_raw)))

## ---- 6. tree table in makeObsData()'s schema -------------------------------
trees_raw <- trees_raw[STATUSCD %in% c(1L, 2L)]           # no cut (3) or unsampled records
trees_raw[, treeName := paste(plotID, SUBP, TREE, sep = "_")]
tree_dt <- trees_raw[, .(siteName = plotID, patchName = subplotID, treeName, year = INVYR, SPCD, species_name = SCIENTIFIC_NAME,
                         dbh = dbh_cm, STATUSCD, PREV_TRE_CN, AGENTCD)]
## v0.2 recruit rule: a living tree is "new" at the first retained record of its
## chronology (the sapling records below 12.7 cm are gone, so a tree tracked from
## the microplot is ingrowth here, not a survivor), "alive" at every later one
setorder(tree_dt, treeName, year)
tree_dt[, first_retained := year == min(year), by = treeName]
tree_dt[, status := fifelse(STATUSCD == 1L & first_retained, "new",
                   fifelse(STATUSCD == 1L, "alive",
                   fifelse(STATUSCD == 2L, "dead", NA_character_)))]
tree_dt[, first_retained := NULL]
tree_dt[, living := status %chin% c("alive", "new")]
tree_dt[, status_before := shift(status, 1L, type = "lag"), by = treeName]
tree_dt[, mort_cause := fifelse(STATUSCD == 2L & AGENTCD == 80L, "land use",
                       fifelse(STATUSCD == 2L & AGENTCD == 30L, "fire", fifelse(STATUSCD == 2L, "natural", NA_character_)))]
land_use_sites <- unique(tree_dt[mort_cause == "land use"]$siteName)
if (length(land_use_sites)) {
  message(sprintf("[filter] %d plots with land-use mortality dropped", length(land_use_sites)))
  tree_dt <- tree_dt[!(siteName %in% land_use_sites)]
}
tree_dt[, mort := status == "dead" & status_before %chin% c("alive", "new")]
tree_dt[, c("STATUSCD", "PREV_TRE_CN", "AGENTCD", "SPCD", "status_before", "mort_cause") := NULL]
message(sprintf("[trees] %d records on %d plots, %d subplots, %d inventory years",
                nrow(tree_dt), uniqueN(tree_dt$siteName), uniqueN(tree_dt[, .(siteName, patchName)]), uniqueN(tree_dt$year)))

## ---- 7. static environment per plot (drops plots outside a WorldClim cell) --
env_static <- unique(env_in_dt[plotID %in% unique(tree_dt$siteName), c("plotID", ENV_COLS), with = FALSE])
setnames(env_static, "plotID", "siteName")
bad_env <- env_static[!stats::complete.cases(env_static[, ENV_COLS, with = FALSE]), siteName]
if (length(bad_env)) {
  message(sprintf("[env] %d plots with a missing static covariate dropped", length(bad_env)))
  env_static <- env_static[!siteName %in% bad_env]; tree_dt <- tree_dt[!siteName %in% bad_env]
}
site_years <- unique(tree_dt[, .(siteName, year)])
site_years <- site_years[siteName %in% env_static$siteName]
tree_dt    <- tree_dt[siteName %in% env_static$siteName]
env_dt <- merge(site_years, env_static, by = "siteName", all.x = TRUE)
setcolorder(env_dt, c("siteName", "year", ENV_COLS))

## ---- 8. completeness flag and makeObsData() ---------------------------------
FINN:::add_site_completeness_flags(env_dt = env_dt, tree_dt = tree_dt)
message(sprintf("[completeness] complete plots: %d", uniqueN(tree_dt[complete == TRUE]$siteName)))
obs_list <- makeObsData(tree_dt = tree_dt[complete == TRUE], plotsize = PARAMS$patch_size_ha, aggregate_by_site = FALSE,
                        Npatches = PARAMS$n_patches, minNyears = PARAMS$min_n_years, fix_period_length = PARAMS$fix_period_length,
                        dbh_growth_thresh = PARAMS$dbh_growth_thresh, Nspecies = NULL, NspeciesQuantile = NULL)
obs_dt_raw  <- obs_list$obs_dt
tree_dt_obs <- obs_list$tree_dt[siteName %in% unique(obs_dt_raw$siteName)]
message(sprintf("[makeObsData] %d observation rows on %d plots", nrow(obs_dt_raw), uniqueN(obs_dt_raw$siteName)))

## ---- 9. integer ids ----------------------------------------------------------
inputs <- resolveSiteIDs(tree_dt = as.data.table(tree_dt_obs), obs_dt = as.data.table(obs_dt_raw),
                         env_dt = as.data.table(env_dt), createInitCohorts = FALSE)
obs_dt      <- inputs$obs_dt[year != 0]
obs_dt_p    <- inputs$obs_dt_patches[year != 0]
env_dt      <- inputs$env_dt
siteID_dt   <- inputs$siteID_dt
species_out <- inputs$species_dt[species_name != "other" | species %in% obs_dt$species]
init_trees  <- inputs$tree_dt[year == 0]
init_trees[living == FALSE | is.na(living), trees := 0]
if (!"trees" %in% names(init_trees)) init_trees[, trees := 1L]

## ---- 9b. annual time axis ----------------------------------------------------
cal <- unique(siteID_dt[, .(siteID, period = year, OrigYear)])[order(siteID, period)]
cal[, `:=`(rel_year = OrigYear - OrigYear[1L], period_length = c(NA_integer_, diff(OrigYear))), by = siteID]
cal[, obs_year := period - 1L]
span <- cal[, .(span = max(rel_year), n_inv = .N), by = siteID]
drop_sites <- span[span > PARAMS$max_span_years, siteID]
if (length(drop_sites)) {
  message(sprintf("[time] %d plots with a span > %d years dropped", length(drop_sites), PARAMS$max_span_years))
  obs_dt <- obs_dt[!siteID %in% drop_sites]; obs_dt_p <- obs_dt_p[!siteID %in% drop_sites]; env_dt <- env_dt[!siteID %in% drop_sites]
  init_trees <- init_trees[!siteID %in% drop_sites]; siteID_dt <- siteID_dt[!siteID %in% drop_sites]
  cal <- cal[!siteID %in% drop_sites]; span <- span[!siteID %in% drop_sites]
}
Tmax <- max(cal$rel_year)
key  <- cal[, .(siteID, year = obs_year, rel_year, period_length)]
annualise <- function(dt) {
  dt <- merge(dt, key, by = c("siteID", "year"))
  dt[, `:=`(growth = (1 + growth)^(1 / period_length) - 1, mort = 1 - (1 - mort)^(1 / period_length))]
  dt[, year := rel_year][, rel_year := NULL]
  dt
}
obs_dt   <- annualise(obs_dt)
obs_dt_p <- annualise(obs_dt_p)
ann <- fread(file.path(INT, "annual_climate_by_plot.csv"))[, c("plotID", "year", PARAMS$annual_vars), with = FALSE]
setnames(ann, c("plotID", "year"), c("siteName", "cal_year"))
xy_static <- unique(env_static[, c("siteName", PARAMS$xy_vars), with = FALSE])
site_y0   <- merge(cal[period == 1L, .(siteID, y0 = OrigYear)], unique(siteID_dt[, .(siteID, siteName)]), by = "siteID")
env_dt <- site_y0[, .(year = seq_len(Tmax), cal_year = y0 + seq_len(Tmax)), by = .(siteID, siteName)]
env_dt <- merge(env_dt, ann, by = c("siteName", "cal_year"), all.x = TRUE)
setorder(env_dt, siteID, year)
for (v in PARAMS$annual_vars) env_dt[, (v) := nafill(get(v), type = "locf"), by = siteID]   # fill forward past the last climate year
no_clim <- env_dt[, .(bad = anyNA(.SD)), .SDcols = PARAMS$annual_vars, by = siteID][bad == TRUE, siteID]
if (length(no_clim)) {
  message(sprintf("[time] %d plots without yearly climate dropped", length(no_clim)))
  obs_dt <- obs_dt[!siteID %in% no_clim]; obs_dt_p <- obs_dt_p[!siteID %in% no_clim]; env_dt <- env_dt[!siteID %in% no_clim]
  init_trees <- init_trees[!siteID %in% no_clim]; siteID_dt <- siteID_dt[!siteID %in% no_clim]
  cal <- cal[!siteID %in% no_clim]; span <- span[!siteID %in% no_clim]
}
env_dt <- merge(env_dt, xy_static, by = "siteName", all.x = TRUE)
ENV_COLS <- c(PARAMS$annual_vars, PARAMS$xy_vars)
env_dt <- env_dt[, c("siteID", "year", ENV_COLS), with = FALSE]
stopifnot(nrow(env_dt) == uniqueN(obs_dt$siteID) * Tmax, !anyNA(env_dt[, ENV_COLS, with = FALSE]))
message(sprintf("[time] %d plots, horizon %d years; inventories per plot %s; intervals (years:n) %s",
                uniqueN(obs_dt$siteID), Tmax, paste(names(table(span$n_inv)), table(span$n_inv), sep = ":", collapse = " "),
                paste(names(table(cal$period_length)), table(cal$period_length), sep = ":", collapse = " ")))

## ---- 9c. size classes, quantiles and the tree-level table ------------------
trees_all <- as.data.table(inputs$tree_dt)[siteID %in% unique(obs_dt$siteID)]
key_year  <- cal[, .(siteID, year = obs_year, rel_year, period_length)]
sp_code   <- unique(species_out[, .(species_name, species)])
brk <- PARAMS$size_class_breaks; K <- length(brk) + 1L; CLS <- paste0("n_class_", seq_len(K))
cls_dt <- trees_all[living == TRUE & year > 0 & !is.na(dbh), .(siteID, year, species, cls = findInterval(dbh, brk) + 1L)]
cls_dt <- cls_dt[, .(n = .N), by = .(siteID, year, species, cls)]
cls_dt <- merge(cls_dt, key_year, by = c("siteID", "year"))[, year := rel_year]
cls_w  <- dcast(cls_dt, siteID + year + species ~ cls, value.var = "n", fill = 0L)
for (k in seq_len(K)) if (!as.character(k) %chin% names(cls_w)) cls_w[[as.character(k)]] <- 0L
setnames(cls_w, as.character(seq_len(K)), CLS)
obs_dt <- merge(obs_dt, cls_w[, c("siteID", "year", "species", CLS), with = FALSE], by = c("siteID", "year", "species"), all.x = TRUE)
for (cc in CLS) obs_dt[is.na(get(cc)), (cc) := 0L]
wq <- function(x, p) { x <- sort(x[is.finite(x)]); if (!length(x)) return(rep(NA_real_, length(p))); x[pmax(1L, ceiling(p * length(x)))] }
QCOLS <- sprintf("dbh_q%g", 100 * PARAMS$dbh_probs)
quant_dt <- trees_all[living == TRUE & year > 0 & !is.na(dbh), as.list(wq(dbh, PARAMS$dbh_probs)), by = .(siteID, year, species)]
setnames(quant_dt, setdiff(names(quant_dt), c("siteID", "year", "species")), QCOLS)
quant_dt <- merge(quant_dt, key_year, by = c("siteID", "year"))[, year := rel_year]
obs_dt <- merge(obs_dt, quant_dt[, c("siteID", "year", "species", QCOLS), with = FALSE], by = c("siteID", "year", "species"), all.x = TRUE)

tree_obs_dt <- trees_all[living_before == TRUE & year > 0 & !is.na(dbh_before),
                         .(siteID, year, patchID, species_name = species_before, dbh = dbh_before, rel_growth, died = as.integer(status == "dead"))]
tree_obs_dt <- merge(tree_obs_dt, sp_code, by = "species_name")
tree_obs_dt <- merge(tree_obs_dt, key_year, by = c("siteID", "year"))
tree_obs_dt[, `:=`(growth = (1 + rel_growth)^(1 / period_length) - 1, year = rel_year, trees = 1L)]
tree_obs_dt <- tree_obs_dt[, .(siteID, year, patchID, species, dbh, growth, died, period_length, trees)]
setorder(tree_obs_dt, siteID, year, patchID)
chk_m <- merge(tree_obs_dt[, .(n_at_risk_tree = .N, n_died_tree = sum(died)), by = .(siteID, year, species)],
               obs_dt[, .(siteID, year, species, n_at_risk, n_died)], by = c("siteID", "year", "species"), all.x = TRUE)
stopifnot(chk_m[, all(n_at_risk_tree >= n_at_risk & n_died_tree >= n_died)],
          chk_m[, sum(n_at_risk_tree - n_at_risk)] <= 0.005 * chk_m[, sum(n_at_risk)],
          chk_m[, sum(n_died_tree - n_died)] <= 0.005 * chk_m[, sum(n_died)])
message(sprintf("[tree obs] %d tree intervals, %d with growth, %d deaths", nrow(tree_obs_dt), sum(is.finite(tree_obs_dt$growth)), sum(tree_obs_dt$died)))

## ---- 10. z-scaled copy of the environment -----------------------------------
env_unscaled_dt <- copy(env_dt)
env_scales_dt <- data.table()
for (v in ENV_COLS) {
  mu <- mean(env_dt[[v]], na.rm = TRUE); sdv <- sd(env_dt[[v]], na.rm = TRUE)
  env_scales_dt <- rbind(env_scales_dt, data.table(variable = v, mean = mu, sd = sdv))
  env_dt[[v]] <- (env_dt[[v]] - mu) / sdv
}
env_keep <- c("siteID", "year", ENV_COLS)
env_dt <- env_dt[, ..env_keep]; env_unscaled_dt <- env_unscaled_dt[, ..env_keep]

## ---- 11. spatial blocks and the sets ---------------------------------------
site_xy <- unique(env_unscaled_dt[, .(siteID, lon = x, lat = y)])
site_xy[, block := paste(floor(lon / PARAMS$block_deg), floor(lat / PARAMS$block_deg))]
blocks_all <- sort(unique(site_xy$block))
set.seed(PARAMS$block_seed)
blocks_holdout <- sample(blocks_all, round(PARAMS$holdout_block_share * length(blocks_all)))
site_xy[, block_role := fifelse(block %in% blocks_holdout, "holdout", "train")]
message(sprintf("[blocks] %d blocks of %g degrees; %d held out (%d plots) / %d retained (%d plots)", length(blocks_all), PARAMS$block_deg,
                length(blocks_holdout), site_xy[block_role == "holdout", .N], length(blocks_all) - length(blocks_holdout), site_xy[block_role == "train", .N]))
all_sites <- sort(unique(obs_dt$siteID))
chosen <- list()
for (set_name in names(PARAMS$sets)) {
  spec <- PARAMS$sets[[set_name]]
  pool <- all_sites
  if (!is.null(spec$blocks)) pool <- intersect(pool, site_xy[block_role == spec$blocks, siteID])
  ids <- if (is.null(spec$n)) pool else { set.seed(spec$seed); sort(sample(pool, min(spec$n, length(pool)))) }
  chosen[[set_name]] <- ids
  map <- data.table(siteID = ids, new_id = seq_along(ids))
  renumber <- function(dt) { dt <- merge(dt[siteID %in% ids], map, by = "siteID"); dt[, siteID := new_id][, new_id := NULL]; setcolorder(dt, "siteID"); dt[] }
  out_dir <- file.path(OUT_PREP, set_name); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(renumber(obs_dt), file.path(out_dir, "obs_dt.csv"))
  fwrite(renumber(obs_dt_p), file.path(out_dir, "obs_dt_patches.csv"))
  fwrite(renumber(tree_obs_dt), file.path(out_dir, "tree_obs_dt.csv"))
  fwrite(renumber(env_dt), file.path(out_dir, "env_dt.csv"))
  fwrite(renumber(env_unscaled_dt), file.path(out_dir, "env_unscaled_dt.csv"))
  fwrite(env_scales_dt, file.path(out_dir, "env_scales.csv"))
  fwrite(renumber(init_trees), file.path(out_dir, "init_cohort_dt.csv"))
  fwrite(species_out, file.path(out_dir, "species_dt.csv"))
  fwrite(renumber(siteID_dt), file.path(out_dir, "siteID_dt.csv"))
  fwrite(renumber(cal), file.path(out_dir, "site_calendar.csv"))
  message(sprintf("[set] %-16s %5d plots -> %s", set_name, length(ids), out_dir))
}
tr <- site_xy[siteID %in% chosen$train1500]; ho <- site_xy[siteID %in% chosen$holdout_blocked]
d_km <- vapply(seq_len(nrow(ho)), function(i) {
  dx <- (tr$lon - ho$lon[i]) * 111.32 * cos(ho$lat[i] * pi / 180); dy <- (tr$lat - ho$lat[i]) * 110.57
  min(sqrt(dx^2 + dy^2))
}, numeric(1))
fwrite(data.table(siteID = ho$siteID, lon = ho$lon, lat = ho$lat, nearest_train_km = d_km), file.path(OUT_PREP, "holdout_separation.csv"))
fwrite(site_xy, file.path(OUT_PREP, "site_blocks.csv"))
writeLines(c(
  sprintf("# Prepared data: %s", PARAMS$prep_version), "",
  sprintf("Written %s by `calibration/conus/08_prepare_annual_v5.R` (FINN %s), the port of FINN-fia scripts/07_prepare_finn_inputs.R at commit 58fcf3e.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), as.character(packageVersion("FINN"))), "",
  "Annual time mode with per-plot intervals; yearly TerraClimate tavg_c / prec_mm / cwd_z + x, y; top 20 species of the step-7 ranking,",
  "so the 99 % basal-area rule keeps only plots almost entirely made of modelled species; macroplot thinning; patch = one subplot",
  "(0.0168 ha); recruitment = ingrowth into 12.7 cm; growth / mortality annualised; env years start at 1; train and holdout spatially",
  "blocked (2-degree blocks, 20 % held out); plus size-class stem counts, diameter quantiles and a tree-level observation table.", "",
  "## Sets", "",
  sprintf("- %s: %d plots", names(chosen), lengths(chosen)), "",
  sprintf("Spatial blocking: %g degree blocks, %d of %d blocks held out (seed %d).", PARAMS$block_deg, length(blocks_holdout), length(blocks_all), PARAMS$block_seed), "",
  sprintf("Distance from each holdout plot to the nearest training plot: median %.0f km, minimum %.0f km, 10 %% below %.0f km.",
          median(d_km), min(d_km), quantile(d_km, 0.1)), "",
  sprintf("Species: %d. Env variables: %s.", uniqueN(species_out$species), paste(ENV_COLS, collapse = ", ")), "",
  "## PARAMS", "", "```r", capture.output(str(PARAMS, give.attr = FALSE)), "```"
), file.path(OUT_PREP, "PREPARATION.md"))
message(sprintf("done: %s (%d species, obs years %d..%d, env years 1..%d)", OUT_PREP, uniqueN(species_out$species), min(obs_dt$year), max(obs_dt$year), Tmax))
