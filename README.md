# FINN-model-FIA

Fitted [FINN](https://github.com/FINNverse/FINN) forest models calibrated to the US
Forest Inventory and Analysis (FIA) network, packaged as a versioned simulation tool
(R package `FINNmodelFIA`). The model card is in the
[FINN-zoo](https://finnverse.github.io/FINN-zoo/models.html).

**Current release: v0.1.0** (see `NEWS.md`). The fitted parameters are release
assets of this repository; `load_model()` downloads them once and checks their sha256.

Two variants:

| variant | region | species | time step | climate inputs |
|---|---|---|---|---|
| `oregon` | Oregon, 708 FIA plots | 27 + other | 10 years | `temp`, `prec` (WorldClim 2.1 normals) |
| `conus` | conterminous US, 2500 FIA plots | 20 | 1 year | `tavg_c`, `prec_mm`, `cwd_z` (yearly TerraClimate) |

```r
remotes::install_github("FINNverse/FINN@0848b60")           # FINN 0.2.0, the version the models were fitted with
remotes::install_github("FINNverse/FINN-model-FIA@v0.1.0")
library(FINNmodelFIA)
torch::torch_set_num_threads(1L)   # pinned: the fits and the validation ran single-threaded
m = load_model("conus")            # downloads the released parameters once, checks sha256, rebuilds with the installed FINN
species_table("conus")
model_info()
vignette("quickstart", package = "FINNmodelFIA")
```

Simulations are run with one torch thread and an explicit `batchsize` (250 sites per
`simulateForest()` batch, the value pinned in `model_info()`; `simulate_model()` uses it
by default). A different batchsize changes the random draws but not the model.

This repository is a tool, not the code of a manuscript: it is versioned for users
(see NEWS.md) and may differ from the models in our papers. The calibration pipeline
that produced each release is in `calibration/` (data preparation, fit, validation,
provenance of every raw layer, and the summary of the prepared data a rerun is diffed
against).

License: GPL (>= 3) for code and fitted parameters; the calibration data keep the
license of their source (FIA: US Government work, public domain / CC0 1.0; TerraClimate:
public domain), see `calibration/PROVENANCE.md`. Climate normals: WorldClim 2.1 (Fick &
Hijmans 2017, doi:10.1002/joc.5086), free for academic and other non-commercial use; the
rasters are not redistributed, only values extracted at the plot locations are included.
The raw data are not redistributed here; `calibration/` downloads them from the source.
