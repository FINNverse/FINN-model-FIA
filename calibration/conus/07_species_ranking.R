#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/conus/07_species_ranking.R   (variant "conus", step 7 of 9)
##
## WHAT   Ranks the species of the cohort by the number of plots they are
##        present on, counted on FINN's own criteria (alive, dbh >= 12.7 cm, on
##        plots with four consistently sampled subplots, macroplot records
##        thinned), sweeps top-k species sets under the 99 % basal-area purity
##        rule of step 8, and writes the ranked list up to the "knee" (the k
##        that maximises observations per parameter while every species keeps
##        >= 20 plots). Step 8 takes the top 20 of this list.
## WHY    Port of FINN-fia `scripts/07b_species_scaling.R` (annual-v5 read its
##        `selected_species.csv` of 2026-09-08 and kept the top 20 by rank).
##        The ranking is what step 8 consumes; the sweep table is kept as the
##        record of why 20 species.
## INPUT  <DATA_DIR>/intermediate/conus/{calibration_set_plots.csv, calibration_set_trees.csv}
## OUTPUT <DATA_DIR>/intermediate/conus/selected_species.csv      rank, SPCD, SCIENTIFIC_NAME, n_plots_present (ranks 1..knee)
##        <DATA_DIR>/intermediate/conus/species_scaling_sweep.csv  one row per candidate k
## RUN    Rscript calibration/conus/07_species_ranking.R   (from the repository root, ~2 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))

DATA_DIR <- {
  cands <- c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))
  d <- cands[dir.exists(cands)][1]
  if (is.na(d)) stop("FINN-zoo-data/fia not found (expected two levels above the repository, see calibration/README.md)")
  normalizePath(d)
}
INT <- file.path(DATA_DIR, "intermediate", "conus")

PARAMS <- list(
  min_dbh_cm            = 12.7,    # tally threshold of the FIA subplot; saplings dropped as in step 8
  min_selected_ba_share = 0.99,    # the purity rule of step 8
  k_step                = 5L,      # candidate set sizes: every 5th k plus the total
  min_species_plots     = 20L,     # adequacy floor: every selected species on >= 20 kept plots
  params_per_species    = 8L,      # nominal parameters per species (scales obs/param; the knee does not depend on it)
  macroplot_dbh_keep_one_in = 6L   # Pacific-Northwest macroplot records (>= 61 cm) represent 1/6 of the subplot density
)

## macroplot thinning, identical to step 8: keep one in six macroplot records per
## physical tree, by a deterministic hash of the tree id
thin_macroplot_trees <- function(dt, keep_one_in) {
  is_macr <- dt$geom == "macr"
  if (!any(is_macr, na.rm = TRUE)) return(dt)
  key <- paste(dt$plotID, dt$SUBP, dt$TREE, sep = "_")[which(is_macr)]
  h <- vapply(unique(key), function(s) { v <- utf8ToInt(s); sum(v * seq_along(v)) %% keep_one_in }, numeric(1))
  keep <- h[key] == 0
  drop_idx <- which(is_macr)[!keep]
  message(sprintf("[macroplot] %d macroplot records on %d plots; keeping 1 in %d per tree -> dropped %d",
                  sum(is_macr, na.rm = TRUE), uniqueN(dt$plotID[which(is_macr)]), keep_one_in, length(drop_idx)))
  if (length(drop_idx)) dt <- dt[-drop_idx]
  dt
}

plots_dt <- fread(file.path(INT, "calibration_set_plots.csv"))
plots_dt <- plots_dt[n_subp_sampled_first == 4L & subplots_consistent_over_time == TRUE]
pilot_plotIDs <- plots_dt$plotID
message(sprintf("[plots] %d plots after the 4-subplot filter", length(pilot_plotIDs)))
trees <- fread(file.path(INT, "calibration_set_trees.csv"),
               select = c("plotID", "SUBP", "TREE", "geom", "INVYR", "SPCD", "SCIENTIFIC_NAME", "STATUSCD", "dbh_cm", "ba_m2"))
trees <- thin_macroplot_trees(trees, PARAMS$macroplot_dbh_keep_one_in)
trees <- trees[plotID %in% pilot_plotIDs & STATUSCD == 1L & !is.na(dbh_cm) & dbh_cm >= PARAMS$min_dbh_cm]
message(sprintf("[trees] %d alive records (dbh >= %.1f cm) on %d plots", nrow(trees), PARAMS$min_dbh_cm, uniqueN(trees$plotID)))

sp_rank <- trees[, .(n_plots_present = uniqueN(plotID)), by = .(SPCD, SCIENTIFIC_NAME)][order(-n_plots_present)]
sp_rank[, rank := .I]
M <- nrow(sp_rank)
trees <- merge(trees, sp_rank[, .(SPCD, rank)], by = "SPCD")

first_yr <- trees[, .(first_yr = min(INVYR)), by = plotID]
ba_first <- merge(trees[, .(plotID, INVYR, rank, ba_m2)], first_yr, by = "plotID")[INVYR == first_yr]
ba_pr <- ba_first[, .(ba = sum(ba_m2, na.rm = TRUE)), by = .(plotID, rank)]
setorder(ba_pr, plotID, rank)
ba_pr[, ba_tot := sum(ba), by = plotID]
ba_pr[, cum_share := cumsum(ba) / ba_tot, by = plotID]
required_k <- ba_pr[cum_share >= PARAMS$min_selected_ba_share, .(required_k = min(rank)), by = plotID]
required_k <- merge(data.table(plotID = pilot_plotIDs), required_k, by = "plotID", all.x = TRUE)
required_k[is.na(required_k), required_k := Inf]
pres <- merge(unique(trees[, .(plotID, rank)]), required_k, by = "plotID")
tree_rank_ct <- merge(trees[, .N, by = .(plotID, rank)], required_k, by = "plotID")

kset <- sort(unique(c(seq(PARAMS$k_step, M, by = PARAMS$k_step), M)))
sweep <- rbindlist(lapply(kset, function(k) {
  n_sites <- required_k[required_k <= k, .N]
  pc <- pres[rank <= k & required_k <= k]
  sp_plots <- pc[, .(nplots = uniqueN(plotID)), by = rank]
  min_species_plots <- if (nrow(sp_plots) < k) 0L else min(sp_plots$nplots)
  data.table(k = k, n_sites = n_sites, total_presence = nrow(pc),
             total_tree_obs = tree_rank_ct[rank <= k & required_k <= k, sum(N)],
             n_present_species = nrow(sp_plots), min_species_plots = min_species_plots,
             med_species_plots = if (nrow(sp_plots)) as.numeric(median(sp_plots$nplots)) else 0,
             n_params = k * PARAMS$params_per_species,
             obs_per_param = nrow(pc) / (k * PARAMS$params_per_species),
             passes_floor = min_species_plots >= PARAMS$min_species_plots)
}))
fwrite(sweep, file.path(INT, "species_scaling_sweep.csv"))
ok <- sweep[passes_floor == TRUE]
knee_k <- if (nrow(ok)) ok[which.max(obs_per_param), k] else sweep[which.max(obs_per_param), k]
message(sprintf("[sweep] %d species present; knee (max obs/param within the floor of %d plots) at k = %d: %d sites",
                M, PARAMS$min_species_plots, knee_k, sweep[k == knee_k, n_sites]))
sel <- sp_rank[rank <= knee_k, .(rank, SPCD, SCIENTIFIC_NAME, n_plots_present)]
fwrite(sel, file.path(INT, "selected_species.csv"))
message(sprintf("[out] selected_species.csv: %d ranked species; top 20: %s", nrow(sel),
                paste(head(sel$SCIENTIFIC_NAME, 20), collapse = ", ")))
