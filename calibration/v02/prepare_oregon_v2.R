#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/prepare_oregon_v2.R   (variant "oregon", data v2)
##
## WHAT   The v0.2 Oregon data (oregon-10yr-v2), built from the same Oregon
##        export as v1 (calibration/oregon/01_prepare.R) under the agreed
##        rare-species rule (V02_BRIEF.md, 05:40 UTC) instead of a lumped
##        `other` class:
##        - the species set S = the 18 species of calibration/v02/species_set_oregon.csv
##          (species_tradeoff.R: of the 27 export species those with >= 10 plots,
##          >= 50 trees at the start census, presence in >= 1 held-out block and
##          a temperature span >= 10 % of the pool's; a coverage choice);
##        - a plot is kept when species outside S hold < 2 % of its basal area
##          at every census that STARTS an interval (steps 0 and 1); the end
##          census (step 2) is never a condition, its out-of-S share is only
##          recorded (siteID_dt.csv: share_out_y0, share_out_y1, share_out_end)
##          so the validation can stratify by it;
##        - trees of species outside S on a kept plot are dropped from states,
##          rates and recruitment (option (a) of the rule; x = 2 % keeps the
##          competition error small).
##        Everything else as v1: macroplot thinning, FINN::makeObsData() per
##        subplot, 1-degree blocks with 6 of 28 held out (seed 21), plots
##        without a living tree at step 0 dropped.
## INPUT  <DATA_DIR>/raw/oregon-finn-data-raw/{full_tree_dt.csv, env_unscaled_dt.csv, oregon_plot_coords.csv},
##        calibration/v02/species_set_oregon.csv
## OUTPUT <DATA_DIR>/prepared/oregon-10yr-v2/{all,train_blocked,holdout_blocked}/ (v1 layout; species_dt has 18 rows;
##        siteID_dt.csv carries the out-of-S shares), PREPARATION.md, site_blocks.csv, holdout_separation.csv,
##        calibration/v02/expected_preparation_oregon_v2.csv
## RUN    Rscript calibration/v02/prepare_oregon_v2.R   (repo root, ~1 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
.libPaths(c("~/R/lib-finn-0848b60", .libPaths()))
suppressPackageStartupMessages({ library(FINN); library(data.table) })
DATA_DIR <- normalizePath(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))[
  dir.exists(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia")))][1])

PARAMS <- list(
  prep_version    = "oregon-10yr-v2",
  raw_dir         = file.path(DATA_DIR, "raw", "oregon-finn-data-raw"),
  species_file    = file.path("calibration", "v02", "species_set_oregon.csv"),   # S, from species_tradeoff.R
  out_share_max   = 0.02,          # species outside S hold < 2 % of the basal area at every interval-start census
  patch_size_ha   = 0.01681, n_patches = 4L,
  macroplot_dbh   = 61, macroplot_keep = 6L,   # Pacific-Northwest macroplot records: 1 in 6 per tree
  min_n_years     = 2L, fix_period = 1L, dbh_growth_thresh = c(-10, 50),
  block_deg       = 1, holdout_block_share = 0.2, block_seed = 21L,
  env_vars        = c("temp", "prec")
)
OUT <- file.path(DATA_DIR, "prepared", PARAMS$prep_version)
if (dir.exists(OUT)) unlink(OUT, recursive = TRUE)        # v2 was never used by a fit; this selection replaces it
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- 1. raw tables and the species set -------------------------------------
trees   <- fread(file.path(PARAMS$raw_dir, "full_tree_dt.csv"))
env_raw <- unique(fread(file.path(PARAMS$raw_dir, "env_unscaled_dt.csv")), by = c("siteID", "year"))
SPECIES <- fread(PARAMS$species_file)[, .(species, species_name)][order(species)]   # named in capitals: `species` is a column of the tree table
stopifnot(identical(SPECIES$species, seq_len(nrow(SPECIES))))
coords <- fread(file.path(PARAMS$raw_dir, "oregon_plot_coords.csv"))[, siteName := paste(COUNTYCD, PLOT, sep = "_")]
sites <- merge(unique(trees[, .(siteName, siteID_raw = siteID)]), coords[, .(siteName, lon, lat)], by = "siteName", all.x = TRUE)
if (anyNA(sites$lat)) stop("plots without coordinates")
message(sprintf("raw: %d tree records, %d plots; S = %d species", nrow(trees), uniqueN(trees$siteName), nrow(SPECIES)))

## ---- 2. macroplot thinning (as v1) -------------------------------------------
tree_key <- trees[, .(macro = any(dbh >= PARAMS$macroplot_dbh, na.rm = TRUE)), by = treeName]
hash_keep <- function(s, k) vapply(s, function(x) { v <- utf8ToInt(x); (sum(v * seq_along(v)) %% k) == 0 }, logical(1))
tree_key[, keep := !macro]; tree_key[macro == TRUE, keep := hash_keep(treeName, PARAMS$macroplot_keep)]
trees <- trees[treeName %in% tree_key[keep == TRUE, treeName]]

## ---- 3. the site rule on the interval-start censuses -------------------------
ba_c <- trees[living %in% TRUE, .(ba = sum(pi * (dbh / 200)^2)), by = .(siteName, year, inS = species_name %in% SPECIES$species_name)]
share <- dcast(ba_c, siteName + year ~ inS, value.var = "ba", fill = 0)
setnames(share, c("FALSE", "TRUE"), c("ba_out", "ba_in"), skip_absent = TRUE)
if (!"ba_out" %in% names(share)) share[, ba_out := 0]
share[, share_out := fifelse(ba_in + ba_out > 0, ba_out / (ba_in + ba_out), 0)]
sh <- dcast(share, siteName ~ year, value.var = "share_out", fill = 0)
setnames(sh, c("0", "1", "2"), c("share_out_y0", "share_out_y1", "share_out_end"))
keep_rule <- sh[share_out_y0 < PARAMS$out_share_max & share_out_y1 < PARAMS$out_share_max, siteName]
message(sprintf("site rule (< %.0f %% out-of-S basal area at steps 0 and 1): %d of %d plots kept; end-census share of the kept plots: < 1 %%: %d, 1-5 %%: %d, > 5 %%: %d",
                100 * PARAMS$out_share_max, length(keep_rule), nrow(sh),
                sh[siteName %in% keep_rule & share_out_end < 0.01, .N], sh[siteName %in% keep_rule & share_out_end >= 0.01 & share_out_end <= 0.05, .N],
                sh[siteName %in% keep_rule & share_out_end > 0.05, .N]))
trees <- trees[siteName %in% keep_rule & species_name %in% SPECIES$species_name]     # option (a): out-of-S trees dropped

## ---- 4. observations ----------------------------------------------------------
tree_dt <- trees[, .(siteName, patchName, treeName, year, species_name, dbh, status, living = living %in% TRUE)]
obs_list <- makeObsData(tree_dt, plotsize = PARAMS$patch_size_ha, aggregate_by_site = TRUE, minNyears = PARAMS$min_n_years,
                        fix_period_length = PARAMS$fix_period, dbh_growth_thresh = PARAMS$dbh_growth_thresh, Npatches = PARAMS$n_patches)
obs <- as.data.table(obs_list$obs_dt)[year > 0]
obs <- merge(obs, SPECIES, by = "species_name", all.x = TRUE)
if (anyNA(obs$species)) stop("species without a code")
kept_sites <- sort(unique(obs$siteName))
grid <- CJ(siteName = kept_sites, year = sort(unique(obs$year)), species = SPECIES$species)
obs <- merge(grid, obs, by = c("siteName", "year", "species"), all.x = TRUE)
obs <- merge(obs[, species_name := NULL], SPECIES, by = "species")
for (cc in c("ba", "trees", "reg")) obs[is.na(get(cc)), (cc) := 0]
for (cc in c("n_at_risk", "n_died", "growth_n")) obs[is.na(get(cc)), (cc) := 0L]
init <- merge(trees[year == 0 & living %in% TRUE & siteName %in% kept_sites, .(siteName, patchName, patchID, treeName, species_name, dbh, trees = 1L)], SPECIES, by = "species_name")
empty0 <- setdiff(kept_sites, unique(init$siteName))
if (length(empty0)) { message(sprintf("%d plots without a living tree at step 0 dropped", length(empty0))); kept_sites <- setdiff(kept_sites, empty0); obs <- obs[siteName %in% kept_sites] }
message(sprintf("makeObsData: %d plots kept (%d of the raw export), %d observation rows", length(kept_sites), uniqueN(sites$siteName), nrow(obs)))

## ---- 5. blocks and the split (as v1) ------------------------------------------
sites <- sites[siteName %in% kept_sites]
sites[, block := paste(floor(lon / PARAMS$block_deg), floor(lat / PARAMS$block_deg))]
blocks <- sites[, .N, by = block][order(block)]
set.seed(PARAMS$block_seed); n_hold <- max(1L, round(PARAMS$holdout_block_share * nrow(blocks))); hold_blocks <- sample(blocks$block, n_hold)
sites[, block_role := fifelse(block %in% hold_blocks, "holdout", "train")]
hav <- function(lon1, lat1, lon2, lat2) { r <- pi / 180; dlat <- (lat2 - lat1) * r; dlon <- (lon2 - lon1) * r
  a <- sin(dlat / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin(dlon / 2)^2; 2 * 6371 * asin(pmin(1, sqrt(a))) }
tr <- sites[block_role == "train"]; ho <- sites[block_role == "holdout"]
ho[, nearest_train_km := vapply(seq_len(.N), function(i) min(hav(lon[i], lat[i], tr$lon, tr$lat)), numeric(1))]
sep <- sprintf("median %.0f km, minimum %.0f km, 10 %% below %.0f km", median(ho$nearest_train_km), min(ho$nearest_train_km), quantile(ho$nearest_train_km, 0.1))
message(sprintf("blocks: %d of %d held out -> %d train, %d holdout plots; separation %s", n_hold, nrow(blocks), nrow(tr), nrow(ho), sep))
fwrite(sites[, .(siteName, lon, lat, block, block_role)], file.path(OUT, "site_blocks.csv"))
fwrite(ho[, .(siteName, lon, lat, nearest_train_km)], file.path(OUT, "holdout_separation.csv"))

## ---- 6. the sets ----------------------------------------------------------------
write_set <- function(set_sites, name) {
  d <- file.path(OUT, name); dir.create(d, showWarnings = FALSE)
  key <- data.table(siteName = sort(set_sites))[, siteID := .I]
  key <- merge(key, sites[, .(siteName, siteID_raw, lon, lat, block, block_role)], by = "siteName")
  key <- merge(key, sh[, .(siteName, share_out_y0, share_out_y1, share_out_end)], by = "siteName")
  o <- merge(obs[siteName %in% set_sites], key[, .(siteName, siteID)], by = "siteName")[
    , .(siteID, year, species, ba, dbh, trees, growth, growth_n, mort, n_at_risk, n_died, reg, species_name)]
  setorder(o, siteID, year, species)
  e <- merge(env_raw[year >= 1], key[, .(siteID_raw, siteID)], by.x = "siteID", by.y = "siteID_raw")
  e <- e[, c("siteID.y", "year", PARAMS$env_vars), with = FALSE]; setnames(e, "siteID.y", "siteID"); setorder(e, siteID, year)
  stopifnot(all(unique(o[, .(siteID, year)])[, paste(siteID, year)] %in% e[, paste(siteID, year)]))
  i <- merge(init[siteName %in% set_sites], key[, .(siteName, siteID)], by = "siteName")[, .(siteID, patchID, species, dbh, trees, treeName, species_name)]
  setorder(i, siteID, patchID, species)
  fwrite(o, file.path(d, "obs_dt.csv")); fwrite(e, file.path(d, "env_dt.csv")); fwrite(i, file.path(d, "init_cohort_dt.csv"))
  fwrite(SPECIES, file.path(d, "species_dt.csv"))
  fwrite(key[order(siteID), .(siteID, siteName, lon, lat, block, block_role, share_out_y0, share_out_y1, share_out_end)], file.path(d, "siteID_dt.csv"))
  message(sprintf("%-16s %4d sites, %6d obs rows, %6d init trees", name, nrow(key), nrow(o), nrow(i)))
  invisible(o)
}
o_all <- write_set(sites$siteName, "all"); o_tr <- write_set(tr$siteName, "train_blocked"); o_ho <- write_set(ho$siteName, "holdout_blocked")

## ---- 7. summary and PREPARATION.md -----------------------------------------------
summarise_set <- function(o, name) {
  tot <- o[, .(set = name, level = "total", species = "all", period = "all", n_sites = uniqueN(siteID), n_rows = .N,
               sum_ba = round(sum(ba), 4), sum_trees = round(sum(trees), 1), sum_reg = round(sum(reg), 2),
               n_at_risk = sum(n_at_risk), n_died = sum(n_died), growth_n = sum(growth_n), mean_growth = round(weighted.mean(growth, growth_n, na.rm = TRUE), 6))]
  sp <- o[, .(level = "species", n_sites = uniqueN(siteID[trees > 0]), n_rows = .N, sum_ba = round(sum(ba), 4), sum_trees = round(sum(trees), 1),
              sum_reg = round(sum(reg), 2), n_at_risk = sum(n_at_risk), n_died = sum(n_died), growth_n = sum(growth_n),
              mean_growth = round(weighted.mean(growth, growth_n, na.rm = TRUE), 6)), by = .(species = species_name, period = as.character(year))][, set := name]
  rbind(tot, sp[, names(tot), with = FALSE])
}
expected <- rbind(summarise_set(o_all, "all"), summarise_set(o_tr, "train_blocked"), summarise_set(o_ho, "holdout_blocked"))[, variant := "oregon_v2"]
fwrite(expected, file.path("calibration", "v02", "expected_preparation_oregon_v2.csv"))
writeLines(c(
  sprintf("# Prepared data: %s", PARAMS$prep_version), "",
  sprintf("Written %s by `calibration/v02/prepare_oregon_v2.R` (FINN %s). Same source, thinning, observation definitions and blocking as oregon-10yr-v1 (its PREPARATION.md); what differs is the species set and the site rule.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), as.character(packageVersion("FINN"))), "",
  "## Species set and site rule (v0.2, agreed rare-species rule)", "",
  sprintf("- S = %d species (`calibration/v02/species_set_oregon.csv`, chosen by `species_tradeoff.R`: of the 27 export species those with >= 10 plots, >= 50 trees at the start census, presence in >= 1 held-out block and a temperature span >= 10 %% of the pool's -- a coverage choice). No lumped class.", nrow(SPECIES)),
  sprintf("- A plot is kept when species outside S hold < %.0f %% of its basal area at every census that starts an interval (steps 0 and 1). The end census is not a condition; its out-of-S share is recorded per plot (`siteID_dt.csv`: share_out_y0, share_out_y1, share_out_end) for the stratified validation.", 100 * PARAMS$out_share_max),
  "- Trees of species outside S on kept plots are dropped from states, rates and recruitment (option (a)).",
  sprintf("- Population claim: %d of the %d plots of the export (%.0f %%) meet the rule.", length(keep_rule), nrow(sh), 100 * length(keep_rule) / nrow(sh)), "",
  "## Sets", "",
  sprintf("- all: %d plots; train_blocked: %d; holdout_blocked: %d (%d of %d 1-degree blocks, seed %d)", uniqueN(o_all$siteID), uniqueN(o_tr$siteID), uniqueN(o_ho$siteID), n_hold, nrow(blocks), PARAMS$block_seed), "",
  sprintf("Distance from each holdout plot to the nearest training plot: %s.", sep), "",
  "## PARAMS", "", "```r", capture.output(str(PARAMS)), "```"
), file.path(OUT, "PREPARATION.md"))
message("done: ", OUT)
