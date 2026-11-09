#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/conus/02_select_cohort.R   (variant "conus", step 2 of 9)
##
## WHAT   Selects the calibration cohort from the plot metadata of step 1: plots
##        whose record begins in 2000 or later, with subplot tree data, at least
##        two inventories on a 5- or 10-year panel, management tier M0 or M1
##        (natural origin, no cutting, no silvicultural mortality), one
##        condition covering >= 95 % of the plot, and whose target species
##        (every species present on >= 100 plots) cover >= 95 % of the plot's
##        basal area. Also writes the target species vector.
## WHY    Port of FINN-fia `scripts/02_select_calibration_set.R` with the
##        parameters of the annual-v5 cohort (min_inventories = 2: the western
##        10-year panel has exactly two visits and FINN handles per-plot interval
##        lengths).
## INPUT  <DATA_DIR>/intermediate/conus/{plot_metadata.csv, plot_species_summary.csv, species_rank.csv}
## OUTPUT <DATA_DIR>/intermediate/conus/calibration_set_plots.csv      one row per cohort plot (full metadata row)
##        <DATA_DIR>/intermediate/conus/calibration_set_species_v.csv  the target species vector V
## RUN    Rscript calibration/conus/02_select_cohort.R    (from the repository root, seconds)
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
  min_inventories       = 2L,                                 # >= 2 visits (the western 10-yr panel has exactly two)
  panel_cycles          = c("5yr", "10yr"),                   # mean remeasurement interval 4-6 or 8-12 years
  management_tiers      = c("M0_pristine", "M1_natural_origin"),   # natural origin, no cut, no silvicultural mortality
  require_post2000      = TRUE,                               # record begins 2000 or later (annualised panel design)
  require_subplot       = TRUE,                               # subplot (>= 12.7 cm) tree data present
  condition_filter      = "cond95",                           # dominant condition covers >= 95 % of the plot
  species_X_pct         = 95L,                                # target species must cover >= 95 % of plot basal area
  min_plots_per_species = 100L,                               # a species enters V if present on >= 100 plots
  states_allow          = NULL                                # all states
)

md  <- fread(file.path(INT, "plot_metadata.csv"))
sp  <- fread(file.path(INT, "plot_species_summary.csv"))
srk <- fread(file.path(INT, "species_rank.csv"))
message(sprintf("[load] %d plots, %d plot-species rows, %d species", nrow(md), nrow(sp), nrow(srk)))

keep <- copy(md)
if (PARAMS$require_post2000) keep <- keep[record_begins_pre2000 == 0L]
if (PARAMS$require_subplot)  keep <- keep[has_subp_dbh5p_0p067ha == 1L]
keep <- keep[n_inventories_raw >= PARAMS$min_inventories]
if (!is.null(PARAMS$panel_cycles))     keep <- keep[panel_cycle %chin% PARAMS$panel_cycles]
if (!is.null(PARAMS$management_tiers)) keep <- keep[management_tier %chin% PARAMS$management_tiers]
keep <- switch(PARAMS$condition_filter, any = keep, cond95 = keep[cond_homogeneous_95 == 1L],
               single_cond = keep[single_condition_plot == 1L], stop("unknown condition_filter"))
if (!is.null(PARAMS$states_allow)) keep <- keep[STATE_ABBR %chin% PARAMS$states_allow]
message(sprintf("[plots] after plot-level filters: %d of %d", nrow(keep), nrow(md)))

V <- srk[n_plots_present >= PARAMS$min_plots_per_species]
message(sprintf("[species] target vector V: %d species present on >= %d plots", nrow(V), PARAMS$min_plots_per_species))

X <- PARAMS$species_X_pct / 100
plot_v_cov <- sp[SPCD %in% V$SPCD & plotID %chin% keep$plotID, .(v_share = sum(ba_share_avg, na.rm = TRUE)), by = plotID]
cohort <- keep[plotID %chin% plot_v_cov[v_share >= X, plotID]]
message(sprintf("[plots] after the %d %% coverage rule: %d plots in %d states", PARAMS$species_X_pct, nrow(cohort), uniqueN(cohort$STATECD)))

fwrite(cohort, file.path(INT, "calibration_set_plots.csv"))
fwrite(V, file.path(INT, "calibration_set_species_v.csv"))
