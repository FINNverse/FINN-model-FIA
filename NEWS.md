# FINNmodelFIA (development)

## v0.1.0 (prepared 2026-11-10, release by the coordinator)

First fitted release. Two variants, both fitted with FINN 0.2.0 (commit 0848b60),
torch 0.15.1, one thread, `weights = "auto"`, `weight_scale = "auto"` (baseline-scaled
losses). Assets per variant: `parameters_<variant>.pt` (named list of parameter tensors,
what `load_model()` rebuilds from) and `model_<variant>.pt` (the whole fitted object);
sha256 of both in `inst/model/weights.json`. Pinned for simulation: one torch thread,
250 sites per `simulateForest()` batch.

* `oregon`: Oregon FIA, 708 plots (590 train, 118 in a spatially blocked holdout of
  1-degree blocks), 27 species + other, 10-year time step, mechanistic processes with a
  linear climate response (`temp`, `prec`; WorldClim 2.1 normals), recruitment with a
  Beverton-Holt cap. Data version `oregon-10yr-v1`. 500 epochs, 18 min.
  Holdout: basal area R2 0.86 (persistence 0.87), stems 0.68 (0.77), growth R2 -0.03
  (Spearman 0.45), mortality AUC 0.60, recruitment R2 -0.41. Bare-ground succession does
  not settle at a plausible basal area (118-415 m2/ha): use for short projections from an
  inventory only.
* `conus`: conterminous US FIA, 1500 training plots and 1000 holdout plots in 2-degree
  blocks, 20 species, annual time step with per-plot remeasurement intervals, climate
  networks (20-20) on yearly TerraClimate `tavg_c`, `prec_mm`, `cwd_z`, recruitment
  limited by conspecific adult density, growth and mortality scored per tree. Data
  version `conus-annual-v1`. 400 epochs, 235 min.
  Holdout: basal area R2 0.82 (persistence 0.82), stems 0.83 (0.83), per-tree growth
  R2 0.27 (Spearman 0.53), per-tree mortality AUC 0.69 (site level 0.68), recruitment
  R2 0.13 (Spearman 0.30). Bare-ground succession converges at 36-56 m2/ha.

The calibration chain (`calibration/`) rebuilds both prepared sets from the public
FIADB download (WorldClim 2.1, TerraClimate), with provenance and checksums of every
raw layer, and reproduces the FINN-fia tables byte for byte.

* Repository created from the FINN model template.
