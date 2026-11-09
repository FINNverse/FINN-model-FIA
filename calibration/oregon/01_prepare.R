#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/oregon/01_prepare.R   (variant "oregon")
##
## WHAT   Builds the FINN input tables of the Oregon FIA model from the
##        tree-level Oregon export that ships with the FINN package sources
##        (FINN/data-raw, written by the FINN-fia analysis pipeline from the
##        public FIADB CSV dump): 775 plots, 4 subplots each, three inventories
##        (2001-03, 2011-13, 2021-24), 27 species plus "other". The model runs
##        on a 10-YEAR TIME STEP: one simulated step is one remeasurement
##        interval, so growth is the relative diameter increment per decade,
##        mortality the fraction of the trees alive at the start of the decade
##        that died, recruitment the ingrowth into the 12.7 cm tally threshold
##        per ha per decade.
## WHY    The bundled FINN::fia_obs_dt is a 200-site, 11-species subsample for
##        the package vignette; the zoo model is calibrated to the whole export.
##        Two things are done differently from the vignette:
##        1. MACROPLOT THINNING. In the Pacific-Northwest FIA design trees
##           >= 61 cm dbh are tallied on the 0.1 ha macroplot around each
##           0.0168 ha subplot, i.e. at 1/6 of the subplot density. The export
##           carries them at full weight, which puts the 90th percentile of
##           basal area at 161 m2/ha on a 0.0168 ha patch. As in FINN-fia
##           (R/fia_macroplot.R) one in six of those trees is kept,
##           deterministically per tree, so every retained record represents
##           the subplot density and the physical subplot area can be used as
##           the FINN patch size.
##        2. SPATIALLY BLOCKED HOLDOUT. 1-degree blocks, 20 % of the blocks held
##           out, training sites only from the retained blocks; the distance of
##           every holdout site to the nearest training site is reported.
## INPUT  <DATA_DIR>/raw/oregon-finn-data-raw/full_tree_dt.csv       tree list, all inventories
##        <DATA_DIR>/raw/oregon-finn-data-raw/env_unscaled_dt.csv    WorldClim 2.1 normals per plot (natural units)
##        <DATA_DIR>/raw/oregon-finn-data-raw/species_dt.csv         species code <-> name
##        FIADB ENTIRE_PLOT.csv (PARAMS$fiadb_plot_csv) for the plot coordinates;
##        read once, the Oregon rows are cached as oregon_plot_coords.csv next
##        to the other raw files so later runs do not need the 450 MB file.
## OUTPUT <DATA_DIR>/prepared/oregon-10yr-v1/{all,train_blocked,holdout_blocked}/
##          obs_dt.csv          site x step x species: ba, dbh, trees, growth, growth_n,
##                              mort, n_at_risk, n_died, reg, species_name
##          env_dt.csv          site x step (years 1..2): temp, prec (+ the other
##                              four WorldClim variables), natural units
##          init_cohort_dt.csv  living trees of the first inventory (step 0)
##          species_dt.csv, siteID_dt.csv (siteID <-> FIA county_plot, lon, lat, block)
##        <DATA_DIR>/prepared/oregon-10yr-v1/PREPARATION.md, site_blocks.csv,
##        holdout_separation.csv
## RUN    Rscript calibration/oregon/01_prepare.R      (from the repository root, ~1 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({
  library(FINN)
  library(data.table)
})

## where the data live (outside git): <parent of the models folder>/FINN-zoo-data/fia
DATA_DIR <- {
  cands <- c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))
  d <- cands[dir.exists(cands)][1]
  if (is.na(d)) stop("FINN-zoo-data/fia not found (expected two levels above the repository, see calibration/README.md)")
  normalizePath(d)
}

PARAMS <- list(
  prep_version    = "oregon-10yr-v1",
  raw_dir         = file.path(DATA_DIR, "raw", "oregon-finn-data-raw"),
  ## the FIADB plot table of the FINN-fia raw dump (public domain, FIA DataMart);
  ## only needed on the first run, the Oregon coordinates are cached
  fiadb_plot_csv  = file.path("..", "..", "FINN-fia", "data", "CSV_FIADB_ENTIRE", "ENTIRE_PLOT.csv"),
  statecd_oregon  = 41L,
  ## FINN patch = one FIA subplot, 24 ft radius = 0.01681 ha (FINN-fia scripts/config.R)
  patch_size_ha   = 0.01681,
  n_patches       = 4L,
  ## Pacific-Northwest design: trees >= 61 cm are macroplot records at 1/6 of
  ## the subplot density (TPA 0.999 vs 6.018); keep one in six per tree
  macroplot_dbh   = 61,
  macroplot_keep  = 6L,
  ## makeObsData: every plot has 4 subplots and 3 inventories, one step apart
  min_n_years     = 2L,
  fix_period      = 1L,
  dbh_growth_thresh = c(-10, 50),     # FINN default: drop sites with an impossible dbh change
  ## spatial blocking: Oregon is ~4 x 8 degrees, so 1-degree blocks give ~30
  ## blocks; 20 % of them held out, drawn with a fixed seed
  block_deg       = 1,
  holdout_block_share = 0.2,
  block_seed      = 21L,
  ## the model's covariates: mean annual temperature and annual precipitation,
  ## the two weakly correlated WorldClim drivers the FINN vignette uses (r = 0.64;
  ## tempmin/tempmax/precseas/precwarmq correlate > 0.86 with them). Only these
  ## two are written, so a user has to supply exactly these two for new sites.
  env_vars        = c("temp", "prec")
)
OUT <- file.path(DATA_DIR, "prepared", PARAMS$prep_version)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- 1. raw tables -----------------------------------------------------------
trees   <- fread(file.path(PARAMS$raw_dir, "full_tree_dt.csv"))
## the export repeats the plot climate once per subplot: one row per plot and step
env_raw <- unique(fread(file.path(PARAMS$raw_dir, "env_unscaled_dt.csv")), by = c("siteID", "year"))
species <- fread(file.path(PARAMS$raw_dir, "species_dt.csv"))
message(sprintf("raw: %d tree records, %d plots, inventories at relative years %s, %d species codes",
                nrow(trees), uniqueN(trees$siteName), paste(sort(unique(trees$year)), collapse = "/"), nrow(species)))

## ---- 2. plot coordinates from FIADB ------------------------------------------
## siteName is "<COUNTYCD>_<PLOT>" within Oregon (STATECD 41); PLOT is unique
## within a county. LAT/LON in FIADB are the public, perturbed coordinates
## (within ~1 mile), which is all that is needed for 1-degree blocks.
coord_file <- file.path(PARAMS$raw_dir, "oregon_plot_coords.csv")
if (!file.exists(coord_file)) {
  message("reading FIADB plot table for the Oregon coordinates: ", PARAMS$fiadb_plot_csv)
  pl <- fread(PARAMS$fiadb_plot_csv, select = c("STATECD", "COUNTYCD", "PLOT", "LAT", "LON"))
  pl <- unique(pl[STATECD == PARAMS$statecd_oregon & is.finite(LAT) & is.finite(LON),
                  .(COUNTYCD, PLOT, LAT, LON)])
  ## one coordinate per plot (identical across inventories; keep the first if not)
  pl <- pl[, .(lat = LAT[1], lon = LON[1]), by = .(COUNTYCD, PLOT)]
  fwrite(pl, coord_file)
}
coords <- fread(coord_file)
coords[, siteName := paste(COUNTYCD, PLOT, sep = "_")]
sites <- unique(trees[, .(siteName, siteID_raw = siteID)])
sites <- merge(sites, coords[, .(siteName, lon, lat)], by = "siteName", all.x = TRUE)
if (anyNA(sites$lat)) stop(sum(is.na(sites$lat)), " plots without coordinates")

## ---- 3. macroplot thinning ---------------------------------------------------
## A tree is a macroplot tree if it reaches the threshold at any inventory; the
## whole history of such a tree is kept or dropped together (otherwise a tree
## would appear and vanish between censuses without a death).
tree_key <- trees[, .(macro = any(dbh >= PARAMS$macroplot_dbh, na.rm = TRUE)), by = treeName]
hash_keep <- function(s, k) {            # deterministic: position-weighted code-point sum
  vapply(s, function(x) { v <- utf8ToInt(x); (sum(v * seq_along(v)) %% k) == 0 }, logical(1))
}
tree_key[, keep := !macro]
tree_key[macro == TRUE, keep := hash_keep(treeName, PARAMS$macroplot_keep)]
n_macro <- sum(tree_key$macro)
trees <- trees[treeName %in% tree_key[keep == TRUE, treeName]]
message(sprintf("macroplot thinning: %d trees >= %g cm, %d kept (1 in %d); %d tree records remain",
                n_macro, PARAMS$macroplot_dbh, tree_key[macro == TRUE, sum(keep)], PARAMS$macroplot_keep, nrow(trees)))

## ---- 4. observations with FINN's own aggregation ----------------------------
## Relative inventory years 0/1/2 are the simulation steps (one step = one
## decade). makeObsData derives growth (dbh/dbh_before - 1 per step), the
## closed-cohort mortality counts, and recruitment per ha from status == "new".
tree_dt <- trees[, .(siteName, patchName, treeName, year, species_name, dbh,
                     status, living = living %in% TRUE)]
obs_list <- makeObsData(tree_dt, plotsize = PARAMS$patch_size_ha, aggregate_by_site = TRUE,
                        minNyears = PARAMS$min_n_years, fix_period_length = PARAMS$fix_period,
                        dbh_growth_thresh = PARAMS$dbh_growth_thresh, Npatches = PARAMS$n_patches)
obs <- as.data.table(obs_list$obs_dt)
obs <- obs[year > 0]                                          # step 0 is the initial state
obs <- merge(obs, species, by = "species_name", all.x = TRUE)
if (anyNA(obs$species)) stop("species without a code: ", paste(unique(obs[is.na(species), species_name]), collapse = ", "))
kept_sites <- sort(unique(obs$siteName))
message(sprintf("makeObsData kept %d of %d plots; %d observation rows at steps %s",
                length(kept_sites), nrow(sites), nrow(obs), paste(sort(unique(obs$year)), collapse = "/")))
## every (site, step) must have a row for every species: fill the zeros
grid <- CJ(siteName = kept_sites, year = sort(unique(obs$year)), species = species$species)
obs <- merge(grid, obs, by = c("siteName", "year", "species"), all.x = TRUE)
obs <- merge(obs[, species_name := NULL], species, by = "species")
for (cc in c("ba", "trees", "reg")) obs[is.na(get(cc)), (cc) := 0]
for (cc in c("n_at_risk", "n_died", "growth_n")) obs[is.na(get(cc)), (cc) := 0L]

## ---- 5. initial cohorts: living trees at step 0 ------------------------------
init <- trees[year == 0 & living %in% TRUE & siteName %in% kept_sites,
              .(siteName, patchName, patchID, treeName, species_name, dbh, trees = 1L)]
init <- merge(init, species, by = "species_name")
## a plot without a living tree at the first inventory has no initial state to
## simulate from (FINN needs at least one cohort per site): dropped
empty0 <- setdiff(kept_sites, unique(init$siteName))
if (length(empty0)) {
  message(sprintf("%d plots without a living tree >= 12.7 cm at the first inventory dropped", length(empty0)))
  kept_sites <- setdiff(kept_sites, empty0)
  obs <- obs[siteName %in% kept_sites]
}

## ---- 6. spatial blocks and the split ----------------------------------------
sites <- sites[siteName %in% kept_sites]
sites[, block := paste(floor(lon / PARAMS$block_deg), floor(lat / PARAMS$block_deg))]
blocks <- sites[, .N, by = block][order(block)]
set.seed(PARAMS$block_seed)
n_hold <- max(1L, round(PARAMS$holdout_block_share * nrow(blocks)))
hold_blocks <- sample(blocks$block, n_hold)
sites[, block_role := fifelse(block %in% hold_blocks, "holdout", "train")]
setorder(sites, block_role, siteName)
## great-circle distance of every holdout site to the nearest training site
hav <- function(lon1, lat1, lon2, lat2) {
  r <- pi / 180; dlat <- (lat2 - lat1) * r; dlon <- (lon2 - lon1) * r
  a <- sin(dlat / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin(dlon / 2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}
tr <- sites[block_role == "train"]; ho <- sites[block_role == "holdout"]
ho[, nearest_train_km := vapply(seq_len(.N), function(i) min(hav(lon[i], lat[i], tr$lon, tr$lat)), numeric(1))]
sep <- sprintf("median %.0f km, minimum %.0f km, 10 %% below %.0f km",
               median(ho$nearest_train_km), min(ho$nearest_train_km), quantile(ho$nearest_train_km, 0.1))
message(sprintf("blocks: %d of %d held out -> %d train sites, %d holdout sites; holdout to nearest training site: %s",
                n_hold, nrow(blocks), nrow(tr), nrow(ho), sep))
fwrite(sites[, .(siteName, lon, lat, block, block_role)], file.path(OUT, "site_blocks.csv"))
fwrite(ho[, .(siteName, lon, lat, nearest_train_km)], file.path(OUT, "holdout_separation.csv"))

## ---- 7. write the sets -------------------------------------------------------
## Each set renumbers its sites 1..n (FINN needs contiguous ids); siteID_dt.csv
## keeps the link to the FIA plot. env rows start at step 1 (no step-0 row):
## simulated step k is driven by env row k and compared with the census at
## step k, exactly as the FINN-fia annual preparation does it.
write_set <- function(set_sites, name) {
  d <- file.path(OUT, name); dir.create(d, showWarnings = FALSE)
  key <- data.table(siteName = sort(set_sites))[, siteID := .I]
  key <- merge(key, sites[, .(siteName, siteID_raw, lon, lat, block, block_role)], by = "siteName")
  o <- merge(obs[siteName %in% set_sites], key[, .(siteName, siteID)], by = "siteName")
  o <- o[, .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, species_name)]
  setorder(o, siteID, year, species)
  e <- merge(env_raw[year >= 1], key[, .(siteID_raw, siteID)], by.x = "siteID", by.y = "siteID_raw")
  e <- e[, c("siteID.y", "year", PARAMS$env_vars), with = FALSE]; setnames(e, "siteID.y", "siteID")
  setorder(e, siteID, year)
  stopifnot(all(unique(o[, .(siteID, year)])[, paste(siteID, year)] %in% e[, paste(siteID, year)]))
  i <- merge(init[siteName %in% set_sites], key[, .(siteName, siteID)], by = "siteName")
  i <- i[, .(siteID, patchID, species, dbh, trees, treeName, species_name)]
  setorder(i, siteID, patchID, species)
  fwrite(o, file.path(d, "obs_dt.csv")); fwrite(e, file.path(d, "env_dt.csv"))
  fwrite(i, file.path(d, "init_cohort_dt.csv")); fwrite(species, file.path(d, "species_dt.csv"))
  fwrite(key[order(siteID), .(siteID, siteName, lon, lat, block, block_role)], file.path(d, "siteID_dt.csv"))
  message(sprintf("%-16s %4d sites, %6d obs rows, %6d init trees -> %s", name, nrow(key), nrow(o), nrow(i), d))
  invisible(o)
}
o_all  <- write_set(sites$siteName, "all")
o_tr   <- write_set(tr$siteName, "train_blocked")
o_ho   <- write_set(ho$siteName, "holdout_blocked")

## ---- 8. expected_preparation.csv: the summary a rerun elsewhere is diffed against
ba_ha <- init[, .(ba = sum(pi * (dbh / 200)^2)), by = .(siteName, patchID)][, ba / PARAMS$patch_size_ha]
sp_tab <- o_all[trees > 0, .(sites_present = uniqueN(siteID)), by = species_name][order(-sites_present)]
summarise_set <- function(o, name) {
  tot <- o[, .(set = name, level = "total", species = "all", period = "all",
               n_sites = uniqueN(siteID), n_rows = .N,
               sum_ba = round(sum(ba), 4), sum_trees = round(sum(trees), 1), sum_reg = round(sum(reg), 2),
               n_at_risk = sum(n_at_risk), n_died = sum(n_died), growth_n = sum(growth_n),
               mean_growth = round(weighted.mean(growth, growth_n, na.rm = TRUE), 6))]
  sp <- o[, .(level = "species", n_sites = uniqueN(siteID[trees > 0]), n_rows = .N,
              sum_ba = round(sum(ba), 4), sum_trees = round(sum(trees), 1), sum_reg = round(sum(reg), 2),
              n_at_risk = sum(n_at_risk), n_died = sum(n_died), growth_n = sum(growth_n),
              mean_growth = round(weighted.mean(growth, growth_n, na.rm = TRUE), 6)),
          by = .(species = species_name, period = as.character(year))][, set := name]
  rbind(tot, sp[, names(tot), with = FALSE])
}
expected <- rbind(summarise_set(o_all, "all"), summarise_set(o_tr, "train_blocked"), summarise_set(o_ho, "holdout_blocked"))
expected[, variant := "oregon"]
init_tab <- init[, .(set = "all", level = "init_trees", species = "all", period = "0", n_sites = uniqueN(siteName),
                     n_rows = .N, sum_ba = round(sum(pi * (dbh / 200)^2), 4), sum_trees = .N, sum_reg = NA_real_,
                     n_at_risk = NA_integer_, n_died = NA_integer_, growth_n = NA_integer_, mean_growth = NA_real_, variant = "oregon")]
expected <- rbind(expected, init_tab)
fwrite(expected, file.path("calibration", "oregon", "expected_preparation.csv"))

## ---- 9. PREPARATION.md -------------------------------------------------------
writeLines(c(
  sprintf("# Prepared data: %s", PARAMS$prep_version), "",
  sprintf("Written %s by `FIA-OR_process/scripts/01_prepare_oregon.R` (FINN %s).",
          format(Sys.time(), "%Y-%m-%d %H:%M"), as.character(packageVersion("FINN"))), "",
  "## Source", "",
  "- Tree list and climate: the Oregon export in `FINN/data-raw/` (`full_tree_dt.csv`, `env_unscaled_dt.csv`,",
  "  `species_dt.csv`; FINN commit 15922a1, 2026-06-24), written by the FINN-fia analysis pipeline",
  "  (`scripts/03_attach_environment.R`, `scripts/07_prepare_finn_inputs.R`) from the FIADB CSV dump",
  "  (USDA Forest Service FIA DataMart, https://apps.fs.usda.gov/fia/datamart/, public domain / CC0; dump of June 2025).",
  "- Plot coordinates: FIADB `ENTIRE_PLOT.csv` (`LAT`, `LON`, the public perturbed coordinates), Oregon rows cached as",
  "  `raw/oregon-finn-data-raw/oregon_plot_coords.csv`.",
  "- Climate: WorldClim 2.1 (1970-2000 normals) at the plot: `temp` (annual mean, C) and `prec` (annual sum, mm),",
  "  constant over the steps, natural units (FINN standardises internally). The export also carries tempmax, tempmin,",
  "  precseas, precwarmq; they correlate > 0.86 with the two kept and are not written.", "",
  "## Definitions", "",
  "- Site = FIA plot, patch = one of its four subplots (0.01681 ha each). One simulation step = one 10-year remeasurement interval;",
  "  steps 1 and 2 are the second and third inventory (2011-13, 2021-24), step 0 (2001-03) is the initial state.",
  sprintf("- Macroplot thinning: %d trees >= %g cm at any inventory, one in %d kept per tree (deterministic hash of the tree id);",
          n_macro, PARAMS$macroplot_dbh, PARAMS$macroplot_keep),
  sprintf("  initial basal area per patch after thinning: median %.0f, 90 %% %.0f, 99 %% %.0f m2/ha.",
          median(ba_ha), quantile(ba_ha, 0.9), quantile(ba_ha, 0.99)),
  "- Observations by `FINN::makeObsData()`: `ba` (m2 per patch, mean over patches), `dbh` (cm, tree-weighted), `trees` (per patch),",
  "  `growth` (relative dbh increment per step, tree-weighted mean, `growth_n` trees behind it), `n_at_risk`/`n_died` (closed cohort,",
  "  `mort = n_died/n_at_risk`), `reg` (recruits into 12.7 cm per ha per step).",
  "- Species: the 27 species of the export plus `other` (code 28), kept as a competitor.", "",
  "## Sets", "",
  sprintf("- all: %d sites", uniqueN(o_all$siteID)),
  sprintf("- train_blocked: %d sites (the %d retained 1-degree blocks)", uniqueN(o_tr$siteID), nrow(blocks) - n_hold),
  sprintf("- holdout_blocked: %d sites (%d of %d blocks, seed %d)", uniqueN(o_ho$siteID), n_hold, nrow(blocks), PARAMS$block_seed), "",
  sprintf("Distance from each holdout site to the nearest training site: %s (`holdout_separation.csv`).", sep), "",
  "## Species (sites where present at any inventory, all sets)", "",
  "| species | sites |", "|---|---|",
  sprintf("| %s | %d |", sp_tab$species_name, sp_tab$sites_present), "",
  "## PARAMS", "", "```r", capture.output(str(PARAMS)), "```"
), file.path(OUT, "PREPARATION.md"))
message("done: ", OUT)
