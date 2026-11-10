# FINNmodelFIA 0.1.0

First fitted release (tag `v0.1.0`). Two variants, both fitted with FINN 0.2.0
(commit 0848b60), torch 0.15.1, one thread, `weights = "auto"`, `weight_scale = "auto"`
(baseline-scaled losses). Assets per variant, attached to the GitHub release:
`parameters_<variant>.pt` (named list of parameter tensors, what `load_model()` rebuilds
from) and `model_<variant>.pt` (the whole fitted object for `torch::torch_load()` under
the same FINN version); sha256 of both in `inst/model/weights.json` (`model_info()`).
Pinned for simulation: one torch thread, 250 sites per `simulateForest()` batch; a
different batchsize changes the random draws but not the model.

* `oregon`: Oregon FIA, 708 plots (590 train, 118 in a spatially blocked holdout of
  1-degree blocks), 27 species + other, 10-year time step, mechanistic processes with a
  linear climate response (`temp`, `prec`; WorldClim 2.1 normals), recruitment with a
  Beverton-Holt cap. Data version `oregon-10yr-v1`. 500 epochs, 18 min.
  Holdout (species x plot): basal area R2 0.86 (persistence of the first inventory
  0.87), stems R2 0.68 (persistence 0.77), species-mean dbh R2 0.50 (0.73); stand totals:
  basal area R2 0.69 (0.68), stems 0.21 (0.53). Growth R2 -0.03 (Spearman 0.45),
  mortality AUC 0.60, recruitment R2 -0.41.
* `conus`: conterminous US FIA, 1500 training plots and 1000 holdout plots in 2-degree
  blocks, 20 species, annual time step with per-plot remeasurement intervals, climate
  networks (20-20) on yearly TerraClimate `tavg_c`, `prec_mm`, `cwd_z`, recruitment
  limited by conspecific adult density, growth and mortality scored per tree. Data
  version `conus-annual-v1`. 400 epochs, 235 min.
  Holdout (species x plot): basal area R2 0.82 (persistence 0.82), stems 0.83 (0.83),
  species-mean dbh 0.55 (0.80); stand totals: basal area R2 0.48 (0.48), stems 0.62
  (0.61). Per-tree growth R2 0.27 (Spearman 0.53), per-tree mortality AUC 0.69 (site
  level 0.68), recruitment R2 0.13 (Spearman 0.30).

Full tables: `calibration/<variant>/metrics.csv`,
`calibration/<variant>/validation/model_vs_persistence.csv`, report in
`calibration/<variant>/report.md`.

Known limitations:

* The states (basal area, stems) are predicted about as well as the persistence
  reference, i.e. most of the holdout skill is the persistence of the start inventory;
  the rates (growth, mortality, recruitment) are weak, recruitment in particular.
* `oregon`: bare-ground succession does not settle at a plausible basal area
  (118-415 m2/ha); use the variant for short projections from an inventory only.
  `conus` converges at 36-56 m2/ha from bare ground.
* FINN's simulated mortality has a floor of 0.01 per model step
  (FINNverse/FINN#51): 1 % per year in `conus`, 1 % per decade in `oregon`. Species
  whose fitted mortality lies below it die at the floor in simulation, which the fit
  did not see.
* The climate response of `conus` is a neural network over three yearly covariates;
  outside the calibration range (the 1958-2023 TerraClimate record at the plots) it
  extrapolates without constraint.

The calibration chain (`calibration/`) rebuilds both prepared sets from the public
FIADB download (WorldClim 2.1, TerraClimate), with provenance and checksums of every
raw layer, and reproduces the FINN-fia tables byte for byte.

Licences: code and fitted parameters GPL (>= 3); FIA data public domain (CC0 1.0),
TerraClimate public domain, WorldClim 2.1 CC BY-SA 4.0 (`calibration/PROVENANCE.md`).

# FINNmodelFIA (development)

* Repository created from the FINN model template.
