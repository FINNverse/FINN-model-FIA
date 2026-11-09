#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/conus/03_extract_trees.R   (variant "conus", step 3 of 9)
##
## WHAT   Extracts the tree records of the cohort plots from the FIADB TREE
##        table, tags each record's sampling geometry (subplot / microplot /
##        macroplot) from its expansion factor, converts units and joins the
##        species names.
## WHY    Port of FINN-fia `scripts/extract_trees.R` (all species, no filter).
##        The 13 GB TREE file is read once here; every later step works on
##        this cohort extract (~450 MB).
## INPUT  <DATA_DIR>/intermediate/conus/calibration_set_plots.csv (step 2);
##        FIADB CSV_FIADB_ENTIRE/ENTIRE_TREE.csv, FIADB_REFERENCE/REF_SPECIES.csv
## OUTPUT <DATA_DIR>/intermediate/conus/calibration_set_trees.csv   one row per tree record
## RUN    Rscript calibration/conus/03_extract_trees.R   (from the repository root, ~5 min, ~8 GB RAM)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))

DATA_DIR <- {
  cands <- c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))
  d <- cands[dir.exists(cands)][1]
  if (is.na(d)) stop("FINN-zoo-data/fia not found (expected two levels above the repository, see calibration/README.md)")
  normalizePath(d)
}
FIADB_DIR <- {
  cands <- c(file.path(DATA_DIR, "raw", "fiadb"), file.path("..", "..", "FINN-fia", "data"))
  d <- cands[dir.exists(file.path(cands, "CSV_FIADB_ENTIRE"))][1]
  if (is.na(d)) stop("FIADB dump not found: expected CSV_FIADB_ENTIRE/ under ", paste(cands, collapse = " or "))
  normalizePath(d)
}
INT <- file.path(DATA_DIR, "intermediate", "conus")

PARAMS <- list(
  tree_cols = c("CN", "PREV_TRE_CN", "STATECD", "UNITCD", "COUNTYCD", "PLOT", "SUBP", "TREE", "INVYR", "SPCD", "STATUSCD",
                "DIA", "HT", "ACTUALHT", "DIAHTCD", "DIACHECK", "AGENTCD", "MORTCD", "STANDING_DEAD_CD", "RECONCILECD",
                "TPA_UNADJ", "CCLCD"),
  tpa_subplot = 6.018046, tpa_microplot = 74.965282, tpa_macroplot = 0.999188   # expansion factors that identify the geometry
)

cohort_ids <- unique(fread(file.path(INT, "calibration_set_plots.csv"), select = "plotID")$plotID)
message(sprintf("[cohort] %d plots", length(cohort_ids)))
t0 <- Sys.time()
tree <- fread(file.path(FIADB_DIR, "CSV_FIADB_ENTIRE", "ENTIRE_TREE.csv"), select = PARAMS$tree_cols, showProgress = FALSE)
message(sprintf("[read] ENTIRE_TREE.csv: %d rows (%.0f s)", nrow(tree), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
tree[, plotID := paste(STATECD, UNITCD, COUNTYCD, PLOT, sep = "_")]
tree[, subplotID := paste(STATECD, UNITCD, COUNTYCD, PLOT, SUBP, sep = "_")]
tree_sel <- tree[plotID %chin% cohort_ids]
rm(tree); invisible(gc(verbose = FALSE))
message(sprintf("[filter] %d tree records on the cohort plots", nrow(tree_sel)))

tree_sel[, geom := fifelse(abs(TPA_UNADJ - PARAMS$tpa_subplot) < 0.01, "subp",
                  fifelse(abs(TPA_UNADJ - PARAMS$tpa_microplot) < 0.5, "micr",
                  fifelse(abs(TPA_UNADJ - PARAMS$tpa_macroplot) < 0.01, "macr", "other")))]
tree_sel[, dbh_cm := DIA * 2.54]
tree_sel[, height_m := HT * 0.3048]
tree_sel[, ba_m2 := pi * (dbh_cm / 200)^2]
tree_sel[, status := fcase(STATUSCD == 1, "alive", STATUSCD == 2, "dead", STATUSCD == 3, "cut", default = NA_character_)]
ref <- fread(file.path(FIADB_DIR, "FIADB_REFERENCE", "REF_SPECIES.csv"),
             select = c("SPCD", "GENUS", "SPECIES", "SCIENTIFIC_NAME", "COMMON_NAME"), showProgress = FALSE)
tree_sel <- merge(tree_sel, ref, by = "SPCD", all.x = TRUE)
out_cols <- intersect(c("plotID", "subplotID", "STATECD", "UNITCD", "COUNTYCD", "PLOT", "SUBP", "TREE", "CN", "PREV_TRE_CN", "INVYR",
                        "SPCD", "SCIENTIFIC_NAME", "COMMON_NAME", "GENUS", "SPECIES", "geom", "TPA_UNADJ", "STATUSCD", "status",
                        "AGENTCD", "MORTCD", "STANDING_DEAD_CD", "RECONCILECD", "DIA", "dbh_cm", "HT", "ACTUALHT", "height_m", "ba_m2",
                        "DIAHTCD", "DIACHECK", "CCLCD"), names(tree_sel))
tree_sel <- tree_sel[, ..out_cols]
setorder(tree_sel, plotID, INVYR, SUBP, TREE)
fwrite(tree_sel, file.path(INT, "calibration_set_trees.csv"))
message(sprintf("[out] calibration_set_trees.csv: %d rows, %d plots, %d species, years %d..%d",
                nrow(tree_sel), uniqueN(tree_sel$plotID), uniqueN(tree_sel$SPCD), min(tree_sel$INVYR), max(tree_sel$INVYR)))
print(tree_sel[, .N, by = geom][order(-N)])
