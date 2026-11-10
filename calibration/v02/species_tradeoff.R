#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
## FINN-model-FIA | calibration/v02/species_tradeoff.R
##
## WHAT   The trade-off behind the v0.2 species set S and the site rule of the
##        agreed rare-species decision (V02_BRIEF.md, 05:40 UTC): a plot is kept
##        when the species outside S hold < x % of its basal area at the START
##        census (never the end census). For both variants and every
##        combination of x in {1, 2, 3, 5} % and candidate |S| (the top-k species
##        by plot presence), the script reports
##          plots kept, kept plots in the held-out blocks, held-out blocks still
##          populated, climate coverage (share of the pool's 2 C x 200 mm climate
##          cells with >= 3 kept plots; the 5-95 % span of temperature and
##          precipitation of the kept plots as a fraction of the pool's),
##          and the identifiability of the species in S on the kept plots:
##          the weakest species' plot count (k), trees at risk at the start
##          census (m), held-out blocks with presence (n) and climate span, and
##          how many species of S fall below k >= 100, m >= 500, n >= 5,
##          span >= 0.5 (the thresholds are a coverage choice, not a measured
##          identifiability limit: 05:50 UTC update).
##        Pools: conus = the 29 847 cohort plots with four consistently sampled
##        subplots (intermediate/conus, steps 1-3, macroplot records thinned 1 in
##        6, trees >= 12.7 cm); Oregon = the 775 plots of the export (same
##        thinning); S candidates: Oregon's 27 export species ranked by presence.
## OUTPUT calibration/v02/species_tradeoff.csv (one row per variant x x x |S|),
##        calibration/v02/species_identifiability_<variant>.csv (per species, for every
##        candidate |S| at x = 2 %), printed tables
## RUN    Rscript calibration/v02/species_tradeoff.R   (repo root, ~2 min)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
suppressPackageStartupMessages(library(data.table))
DATA_DIR <- normalizePath(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia"))[
  dir.exists(c(file.path("..", "..", "FINN-zoo-data", "fia"), file.path("..", "FINN-zoo-data", "fia")))][1])
INT <- file.path(DATA_DIR, "intermediate", "conus")
PARAMS <- list(
  x_pct        = c(1, 2, 3, 5),                    # out-of-S basal-area share at the start census, strict <
  k_conus      = c(10L, 15L, 20L, 25L, 30L, 40L, 50L, 65L),   # candidate |S| by plot presence (conus)
  k_oregon     = c(27L, 24L, 20L, 16L, 12L),       # candidate |S| (Oregon's 27 export species, by presence)
  min_dbh_cm   = 12.7, macro_dbh = 61, macro_keep_one_in = 6L,
  block_deg    = c(conus = 2, oregon = 1), holdout_share = 0.2, block_seed = 21L,   # as the v1 preparations
  climate_bin  = c(tavg_c = 2, prec_mm = 200), min_plots_per_cell = 3L,
  ident        = list(k = 100L, m = 500L, n = 5L, span = 0.5)   # coverage thresholds reported per species
)
thin_macroplot <- function(dt, keep_one_in) {
  is_macr <- dt$geom == "macr"
  if (!any(is_macr, na.rm = TRUE)) return(dt)
  key <- dt$tree_key[which(is_macr)]
  h <- vapply(unique(key), function(s) { v <- utf8ToInt(s); sum(v * seq_along(v)) %% keep_one_in }, numeric(1))
  drop_idx <- which(is_macr)[!(h[key] == 0)]
  if (length(drop_idx)) dt[-drop_idx] else dt
}
blocks_of <- function(sites, deg) {
  sites[, block := paste(floor(lon / deg), floor(lat / deg))]
  b <- sort(unique(sites$block)); set.seed(PARAMS$block_seed)
  hold <- sample(b, round(PARAMS$holdout_share * length(b)))
  sites[, block_role := fifelse(block %in% hold, "holdout", "train")][]
}
## one evaluation: `start` = alive trees at the start census (plotID, species_name, ba, rank, n), `sites` (plotID, lon, lat, tavg_c, prec_mm, block, block_role)
evaluate <- function(start, sites, k, x, variant) {
  inS <- start[, .(ba_tot = sum(ba), ba_out = sum(ba[rank > k])), by = plotID][, share_out := fifelse(ba_tot > 0, ba_out / ba_tot, 1)]
  kept <- inS[share_out < x / 100 & ba_tot > 0, plotID]
  ks <- sites[plotID %in% kept]
  pool_cells <- sites[, uniqueN(paste(floor(tavg_c / PARAMS$climate_bin[["tavg_c"]]), floor(prec_mm / PARAMS$climate_bin[["prec_mm"]])))]
  cells <- ks[, .N, by = .(cell = paste(floor(tavg_c / PARAMS$climate_bin[["tavg_c"]]), floor(prec_mm / PARAMS$climate_bin[["prec_mm"]])))][N >= PARAMS$min_plots_per_cell, .N]
  span <- function(v, pool) diff(quantile(v, c(0.05, 0.95))) / diff(quantile(pool, c(0.05, 0.95)))
  ## identifiability of every species in S on the kept plots
  sp <- merge(start[rank <= k & plotID %in% kept], ks[, .(plotID, block_role, block, tavg_c, prec_mm)], by = "plotID")
  ident <- sp[, .(k_plots = uniqueN(plotID), m_trees = sum(n), n_holdout_blocks = uniqueN(block[block_role == "holdout"]),
                  span_tavg = diff(range(tavg_c)) / diff(range(sites$tavg_c))), by = .(species_name, rank)]
  ident <- merge(data.table(rank = seq_len(k)), ident, by = "rank", all.x = TRUE)
  for (cc in c("k_plots", "m_trees", "n_holdout_blocks")) ident[is.na(get(cc)), (cc) := 0L]
  ident[is.na(span_tavg), span_tavg := 0]
  th <- if (variant == "conus") list(k = 100L, m = 500L, n = 1L, span = 0.2) else list(k = 10L, m = 50L, n = 1L, span = 0.1)   # the coverage choice per pool (CHOICE below)
  fail <- ident[k_plots < th$k | m_trees < th$m | n_holdout_blocks < th$n | span_tavg < th$span, .N]
  list(row = data.table(variant = variant, x_pct = x, n_species = k,
    plots_kept = length(kept), share_of_pool = round(length(kept) / nrow(sites), 3),
    plots_in_holdout_blocks = ks[block_role == "holdout", .N],
    holdout_blocks_populated = ks[block_role == "holdout", uniqueN(block)], holdout_blocks_total = sites[block_role == "holdout", uniqueN(block)],
    climate_cells_covered = round(cells / pool_cells, 3), span_tavg_5_95 = round(span(ks$tavg_c, sites$tavg_c), 3), span_prec_5_95 = round(span(ks$prec_mm, sites$prec_mm), 3),
    ident_min_plots = min(ident$k_plots), ident_min_trees = min(ident$m_trees), ident_min_holdout_blocks = min(ident$n_holdout_blocks), ident_min_span = round(min(ident$span_tavg), 2),
    species_failing = fail), ident = cbind(variant = variant, x_pct = x, n_species = k, ident))
}

## ---- conus ---------------------------------------------------------------------
plots <- fread(file.path(INT, "calibration_set_plots.csv"), select = c("plotID", "lat", "lon", "n_subp_sampled_first", "subplots_consistent_over_time"))
plots <- plots[n_subp_sampled_first == 4L & subplots_consistent_over_time == TRUE]
envp <- fread(file.path(INT, "plot_env.csv"))[, .(plotID, tavg_c, prec_mm)]
sites_c <- merge(plots[, .(plotID, lon, lat)], envp, by = "plotID")[is.finite(tavg_c) & is.finite(prec_mm)]
sites_c <- blocks_of(sites_c, PARAMS$block_deg[["conus"]])
tr <- fread(file.path(INT, "calibration_set_trees.csv"), select = c("plotID", "SUBP", "TREE", "INVYR", "SPCD", "SCIENTIFIC_NAME", "STATUSCD", "geom", "dbh_cm", "ba_m2"))
tr <- tr[plotID %chin% sites_c$plotID & STATUSCD == 1L & !is.na(dbh_cm) & dbh_cm >= PARAMS$min_dbh_cm]
tr[, tree_key := paste(plotID, SUBP, TREE, sep = "_")]
tr <- thin_macroplot(tr, PARAMS$macro_keep_one_in)
tr[, first := min(INVYR), by = plotID]
start_c <- tr[INVYR == first, .(ba = sum(ba_m2), n = .N), by = .(plotID, species_name = SCIENTIFIC_NAME)]
rank_c <- tr[, .(n_plots = uniqueN(plotID)), by = .(species_name = SCIENTIFIC_NAME)][order(-n_plots)][, rank := .I]
start_c <- merge(start_c, rank_c[, .(species_name, rank)], by = "species_name")
message(sprintf("[conus] pool %d plots in %d blocks (%d held out); %d species present at a start census", nrow(sites_c), uniqueN(sites_c$block), sites_c[block_role == "holdout", uniqueN(block)], nrow(rank_c)))
res <- list(); idents <- list()
for (k in PARAMS$k_conus) for (x in PARAMS$x_pct) { e <- evaluate(start_c, sites_c, k, x, "conus"); res[[length(res) + 1]] <- e$row; if (x == 2) idents[[length(idents) + 1]] <- e$ident }
fwrite(rbindlist(idents), file.path("calibration", "v02", "species_identifiability_conus.csv"))
fwrite(rank_c, file.path("calibration", "v02", "species_presence_rank_conus.csv"))

## ---- oregon ----------------------------------------------------------------------
raw <- file.path(DATA_DIR, "raw", "oregon-finn-data-raw")
to <- fread(file.path(raw, "full_tree_dt.csv"))
sp_o <- fread(file.path(raw, "species_dt.csv"))
coords <- fread(file.path(raw, "oregon_plot_coords.csv"))[, siteName := paste(COUNTYCD, PLOT, sep = "_")]
envo <- unique(fread(file.path(raw, "env_unscaled_dt.csv")), by = c("siteID", "year"))[year == 0, .(siteID, tavg_c = temp, prec_mm = prec)]
sites_o <- merge(unique(to[, .(plotID = siteName, siteID)]), coords[, .(plotID = siteName, lon, lat)], by = "plotID")
sites_o <- merge(sites_o, envo, by = "siteID")[, siteID := NULL]
sites_o <- blocks_of(sites_o, PARAMS$block_deg[["oregon"]])
to[, tree_key := treeName][, geom := fifelse(dbh >= PARAMS$macro_dbh, "macr", "subp")]
macro_trees <- to[, .(macro = any(dbh >= PARAMS$macro_dbh, na.rm = TRUE)), by = treeName][macro == TRUE, treeName]
to[, geom := fifelse(treeName %in% macro_trees, "macr", "subp")]
to <- thin_macroplot(to, PARAMS$macro_keep_one_in)
start_o <- to[year == 0 & living %in% TRUE & dbh >= PARAMS$min_dbh_cm, .(ba = sum(pi * (dbh / 200)^2), n = .N), by = .(plotID = siteName, species_name)]
rank_o <- to[living %in% TRUE & species_name != "other", .(n_plots = uniqueN(siteName)), by = species_name][order(-n_plots)][, rank := .I]
rank_o <- rbind(rank_o, data.table(species_name = "other", n_plots = to[species_name == "other", uniqueN(siteName)], rank = nrow(rank_o) + 1L))
start_o <- merge(start_o, rank_o[, .(species_name, rank)], by = "species_name")
message(sprintf("[oregon] pool %d plots in %d blocks (%d held out); `other` ranks last (%d)", nrow(sites_o), uniqueN(sites_o$block), sites_o[block_role == "holdout", uniqueN(block)], nrow(rank_o)))
idents_o <- list()
for (k in PARAMS$k_oregon) for (x in PARAMS$x_pct) { e <- evaluate(start_o, sites_o, k, x, "oregon"); res[[length(res) + 1]] <- e$row; if (x == 2) idents_o[[length(idents_o) + 1]] <- e$ident }
fwrite(rbindlist(idents_o), file.path("calibration", "v02", "species_identifiability_oregon.csv"))
fwrite(rank_o, file.path("calibration", "v02", "species_presence_rank_oregon.csv"))

## ---- the chosen sets: a coverage rule applied to the candidate pools ---------------
## S = the species of the candidate set (conus: top 40 by presence; Oregon: the 27
## export species) that, on the plots kept at x = 2 %, have >= k plots, >= m trees at
## the start census, presence in >= n held-out blocks and a temperature span >= the
## stated fraction of the pool's range. k, m, n, span are a COVERAGE choice (05:50 UTC
## update), scaled to each pool: a species below them cannot get a holdout error of
## its own. The kept plots are then re-evaluated with exactly that S.
CHOICE <- list(conus  = list(candidate = 40L, k = 100L, m = 500L, n = 1L, span = 0.2),
               oregon = list(candidate = 27L, k = 10L,  m = 50L,  n = 1L, span = 0.1))
chosen_rows <- list()
for (v in names(CHOICE)) {
  ch <- CHOICE[[v]]
  idt <- rbindlist(if (v == "conus") idents else idents_o)[n_species == ch$candidate]
  keep_sp <- idt[k_plots >= ch$k & m_trees >= ch$m & n_holdout_blocks >= ch$n & span_tavg >= ch$span][order(rank)]
  dropped <- idt[!species_name %in% keep_sp$species_name][order(rank)]
  message(sprintf("[%s] S = %d of the top %d species; dropped: %s", v, nrow(keep_sp), ch$candidate,
                  paste(sprintf("%s (plots %d, trees %d, holdout blocks %d, span %.2f)", dropped$species_name, dropped$k_plots, dropped$m_trees, dropped$n_holdout_blocks, dropped$span_tavg), collapse = "; ")))
  start <- if (v == "conus") copy(start_c) else copy(start_o); sites <- if (v == "conus") sites_c else sites_o
  ## re-rank so that S occupies ranks 1..|S| and everything else is outside
  start[, rank := fifelse(species_name %in% keep_sp$species_name, match(species_name, keep_sp$species_name), 1000L)]
  for (x in PARAMS$x_pct) {
    e <- evaluate(start, sites, nrow(keep_sp), x, v)
    e$row[, n_species := nrow(keep_sp)][, set := "chosen"]
    chosen_rows[[length(chosen_rows) + 1]] <- e$row
    if (x == 2) fwrite(e$ident[, .(species = rank, species_name, k_plots, m_trees, n_holdout_blocks, span_tavg)],
                       file.path("calibration", "v02", sprintf("species_set_%s.csv", v)))
  }
}
out <- rbind(rbindlist(res)[, set := "candidate"], rbindlist(chosen_rows))
fwrite(out, file.path("calibration", "v02", "species_tradeoff.csv"))
options(width = 220)
print(out[variant == "conus"], nrows = 100); print(out[variant == "oregon"], nrows = 100)
