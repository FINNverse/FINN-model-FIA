# Calibration pipeline

Self-contained R scripts that produced the current release, one per step, each run as
`Rscript <script>` from the repository root with no arguments (house rules:
`r-pipeline-conventions`, `finn-standalone-fit`, `finn-evaluation`). `run_all.R` runs
the whole chain in order, from the public FIADB download to the validated models.

Two model variants are calibrated here:

| variant | region | species | time step | climate | folder |
|---|---|---|---|---|---|
| `oregon` | Oregon, 708 FIA plots (590 train / 118 blocked holdout) | 27 + other | 10 years (one remeasurement interval) | WorldClim normals: `temp`, `prec`, linear per species | `oregon/` |
| `conus` | conterminous US, 2500 FIA plots (1500 train / 1000 blocked holdout) | 20 | 1 year, per-site intervals of 4-12 years | yearly TerraClimate: `tavg_c`, `prec_mm`, `cwd_z`, one 20-20 network per process | `conus/` |

```
calibration/
  run_all.R                   the fixed order of every step
  PROVENANCE.md               every raw layer: source URL, licence, date pulled, sha256
  expected_preparation.csv    summary of the final prepared data of both variants (diff a rerun against it)
  04_package_files.R          writes inst/model/* and inst/extdata/* from the calibration outputs
  hpc/run_fit.sbatch          SLURM wrapper: resources only, calls Rscript <script>
  oregon/
    01_prepare.R              Oregon export + FIADB coordinates -> prepared/oregon-10yr-v1, PREPARATION.md, expected_preparation.csv
    02_fit.R, 03_validation.R
  conus/
    01_fiadb_metadata.R       FIADB dump -> plot metadata, plot x species shares, species ranking
    02_select_cohort.R        the cohort: post-2000, >= 2 inventories, 5/10-yr panel, M0/M1, 95 % condition, 95 % target-species cover
    03_extract_trees.R        the cohort's tree records from ENTIRE_TREE.csv
    04_attach_worldclim.R     WorldClim normals and coordinates per plot
    05_terraclimate_deficit.R TerraClimate def/vpd 1958-latest per plot (OPeNDAP, cached)
    06_terraclimate_annual.R  yearly tavg_c, prec_mm, cwd_z per plot and calendar year
    07_species_ranking.R      species ranked by plot presence under FINN's criteria; the 20 modelled species are its top 20
    08_prepare_annual_v5.R    FINN tables (annual step, per-plot intervals, tree table, 2-degree blocks) -> prepared/annual-v5
    09_prepare_conus.R        reshape for the tool (natural-unit env, no x/y) -> prepared/conus-annual-v1, expected_preparation.csv
    02_fit.R, 03_validation.R
```

Steps `conus/01..08` are the port of the FINN-fia analysis pipeline (scripts 01, 02,
extract_trees, 03, 03c, 03d, 07b, 07 at commit 58fcf3e) with its exact parameters and
seeds; the ported chain reproduces the FINN-fia prepared set `annual-v5` byte for byte
(PROVENANCE.md, "Reproduction check").

## Data location

Raw and prepared data and the fitted `.pt` files are not in git. The scripts expect
them under `../../FINN-zoo-data/fia` relative to the repository root (the layout of the
FINN-zoo working directories; `../FINN-zoo-data/fia` is tried as well):

```
FINN-zoo-data/fia/
  raw/fiadb/{CSV_FIADB_ENTIRE,FIADB_REFERENCE}/   the FIADB dump (PROVENANCE.md 1); a FINN-fia data folder next to the models folder is also found
  raw/worldclim/wc2.1_30s_bio/                    WorldClim 2.1 bio_1, bio_12 (PROVENANCE.md 2); a `worldclim` folder next to the models folder is also found
  raw/oregon-finn-data-raw/                       the Oregon export and the cached coordinates (PROVENANCE.md 4)
  intermediate/conus/                             steps 1-7 of the conus chain, incl. the TerraClimate caches (PROVENANCE.md 3)
  prepared/annual-v5/                             step 8
  prepared/oregon-10yr-v1/, prepared/conus-annual-v1/   the sets the fits read
  fits/<variant>/                                 parameters.pt, model.pt (+ the release-named copies), loss_history.csv, param_history.rds
```

## Software

FINN 0.2.0 at commit 0848b60 (FINNverse/FINN main), installed into its own library
`~/R/lib-finn-0848b60`, which every FINN script puts first on `.libPaths()`. torch 0.15.1
(laptop) / 0.17.0 (cluster); R >= 4.4; data.table, terra (step 4), ncdf4 (steps 5-6),
ggplot2 and gridExtra (validation). All scripts run single-threaded
(`torch_set_num_threads(1)`), so several fits can share a node.

## Weights

`02_fit.R` saves two release assets per variant: `parameters_<variant>.pt`, a **named
list of parameter tensors** that `R/definition.R::build_model()` copies into a freshly
built model (the fit script checks that the rebuilt model reproduces the fitted
simulation exactly), and `model_<variant>.pt`, the whole fitted object for
`torch::torch_load()` under the same FINN version. The standardisation of the climate
inputs (`env_autoscale`) is stored in `inst/model/covariates.csv` (`center`, `scale`).
Loss weighting: `weights = "auto"` with `weight_scale = "auto"` in both fits, i.e. every
loss term is divided by its intercept-only baseline (FINN >= 0.2.0 default); no raw-loss
weights. Pinned for simulation: one torch thread, 250 sites per `simulateForest()` batch
(`weights.json`: `threads`, `sim_batchsize`).
