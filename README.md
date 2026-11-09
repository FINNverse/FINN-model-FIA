# FINN-model-FIA

Fitted [FINN](https://github.com/FINNverse/FINN) forest models calibrated to the US
Forest Inventory and Analysis (FIA) network, packaged as a versioned simulation tool
(R package `FINNmodelFIA`). The model card is in the
[FINN-zoo](https://finnverse.github.io/FINN-zoo/models.html).

**Status: in preparation, no release yet.**

Two variants:

| variant | region | species | time step | climate inputs |
|---|---|---|---|---|
| `oregon` | Oregon, 708 FIA plots | 27 + other | 10 years | `temp`, `prec` (WorldClim 2.1 normals) |
| `conus` | conterminous US, 2500 FIA plots | 20 | 1 year | `tavg_c`, `prec_mm`, `cwd_z` (yearly TerraClimate) |

```r
remotes::install_github("FINNverse/FINN-model-FIA")
library(FINNmodelFIA)
m = load_model("conus")   # downloads the released parameters once, checks sha256, rebuilds with the installed FINN
species_table("conus")
model_info()
vignette("quickstart", package = "FINNmodelFIA")
```

This repository is a tool, not the code of a manuscript: it is versioned for users
(see NEWS.md) and may differ from the models in our papers. The calibration pipeline
that produced each release is in `calibration/` (data preparation, fit, validation,
provenance of every raw layer, and the summary of the prepared data a rerun is diffed
against).

License: GPL (>= 3) for code and fitted parameters; the calibration data keep the
license of their source (FIA: public domain; WorldClim 2.1: CC BY-SA 4.0; TerraClimate:
public domain), see `calibration/PROVENANCE.md`.
