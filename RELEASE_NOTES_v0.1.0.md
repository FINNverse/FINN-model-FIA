First release of FINNmodelFIA: two fitted FINN forest models for US forests, calibrated to the Forest Inventory and Analysis (FIA) network and packaged as an R simulation tool.

Variants: `oregon` (Oregon, 27 species + other, 10-year step, linear climate response on WorldClim 2.1 `temp` and `prec`) and `conus` (conterminous US, 20 species, annual step, climate networks on yearly TerraClimate `tavg_c`, `prec_mm`, `cwd_z`). Both fitted with FINN 0.2.0 (commit 0848b60), torch 0.15.1, one thread, baseline-scaled losses (`weights = "auto"`, `weight_scale = "auto"`).

Install and use:

```r
remotes::install_github("FINNverse/FINN@0848b60")
remotes::install_github("FINNverse/FINN-model-FIA@v0.1.0")
library(FINNmodelFIA)
torch::torch_set_num_threads(1L)
m <- load_model("conus")   # downloads parameters_conus.pt once, checks its sha256, rebuilds the model with the installed FINN
```

Simulations are pinned to one torch thread and `batchsize = 250`; a different batchsize changes the random draws but not the model.

Assets (per variant the name-matched parameter list that `load_model()` rebuilds from, and the whole fitted model object for `torch::torch_load()` under FINN 0848b60):

| file | sha256 |
|---|---|
| parameters_oregon.pt | 0d29d459a7ae10e692550c1ce64d9a94d42d8f03b7a971108ae12c57af32d815 |
| model_oregon.pt | 503b7ecd37ece8cb4144f6db5dad3e9d5c5b0e18af9cf2ef9cd3783b8320092c |
| parameters_conus.pt | d40ce38f94bc307fd17fc172cf12d8c29f5b95415eb08e711d1f1a2838f2a9a9 |
| model_conus.pt | 7c7d41915e64563cf0e2545c94c31c8cf7c0408fc0c27adee84b520f5ae38665 |

Holdout skill (spatially blocked plots the fit never saw; the persistence of the start inventory in brackets): `oregon` basal area R2 0.86 (0.87), stems 0.68 (0.77), mortality AUC 0.60; `conus` basal area R2 0.82 (0.82), stems 0.83 (0.83), per-tree growth R2 0.27, per-tree mortality AUC 0.69. The states are predicted about as well as persistence, the rates are weak; `oregon` is for short projections from an inventory only (its bare-ground succession does not settle). FINN's simulated mortality has a floor of 0.01 per model step (FINNverse/FINN#51), which the fits did not see. Full tables and the limitations are in NEWS.md and `calibration/<variant>/report.md`.

Data: `oregon-10yr-v1` and `conus-annual-v1`, rebuilt end to end from the public FIADB download by `calibration/` (provenance and checksums of every raw layer in `calibration/PROVENANCE.md`). Licences: code and parameters GPL (>= 3); FIA data public domain (CC0 1.0), TerraClimate public domain, WorldClim 2.1 CC BY-SA 4.0.
