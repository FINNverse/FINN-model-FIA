# Calibration pipeline

Self-contained R scripts that produced the current release, one per step, each run as
`Rscript <script>` from the repository root with no arguments (house rules:
`r-pipeline-conventions`, `finn-standalone-fit`, `finn-evaluation`). `run_all.R` runs
the whole chain in order.

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
  oregon/ conus/
    01_prepare.R              raw layers -> prepared/<version>/{train,holdout}, PREPARATION.md, expected_preparation.csv
    02_fit.R                  the production fit; writes parameters.pt (release asset), fit_summary.md, loss_history.csv, env_scaling.csv
    03_validation.R           holdout report: report.md, metrics.csv, validation/*.png
```

## Data location

Raw and prepared data and the fitted `.pt` files are not in git. The scripts expect
them under `../../FINN-zoo-data/fia` relative to the repository root (the layout of the
FINN-zoo working directories; `../FINN-zoo-data/fia` is tried as well):

```
FINN-zoo-data/fia/
  raw/oregon-finn-data-raw/        the Oregon export and the cached coordinates (PROVENANCE.md)
  prepared/annual-v5/              the FINN-fia prepared set the conus variant starts from
  prepared/oregon-10yr-v1/         written by oregon/01_prepare.R
  prepared/conus-annual-v1/        written by conus/01_prepare.R
  fits/<variant>/                  parameters.pt, model.pt, loss_history.csv, param_history.rds
```

## Software

FINN 0.2.0 at commit 0848b60 (FINNverse/FINN main), installed into its own library
`~/R/lib-finn-0848b60`, which every script puts first on `.libPaths()`. torch 0.15.1
(laptop) / 0.17.0 (cluster); R >= 4.4. All scripts run single-threaded
(`torch_set_num_threads(1)`), so several fits can share a node.

## Weights

`02_fit.R` saves the fitted model as a **named list of parameter tensors**
(`parameters.pt`), not as the whole model object; `R/definition.R::build_model()`
recreates the structure with the installed FINN and copies the tensors in, and the fit
script checks that the rebuilt model reproduces the fitted simulation exactly. The
standardisation of the climate inputs (`env_autoscale`) is stored in
`inst/model/covariates.csv` (`center`, `scale`). Loss weighting: `weights = "auto"`
with `weight_scale = "auto"` in both fits, i.e. every loss term is divided by its
intercept-only baseline (FINN >= 0.2.0 default); no raw-loss weights.
