# FINNmodelFIA (development)

## Towards v0.1.0 (in preparation, 2026-11-09)

Two variants, both fitted with FINN 0.2.0 (commit 0848b60), `weights = "auto"`,
`weight_scale = "auto"`:

* `oregon`: Oregon FIA, 708 plots (590 train, 118 in a spatially blocked holdout of
  1-degree blocks), 27 species + other, 10-year time step, mechanistic processes with a
  linear climate response (`temp`, `prec`; WorldClim 2.1 normals), recruitment with a
  Beverton-Holt cap. Data version `oregon-10yr-v1`.
* `conus`: conterminous US FIA, 1500 training plots and 1000 holdout plots in 2-degree
  blocks (FINN-fia `annual-v5`), 20 species, annual time step with per-site
  remeasurement intervals, climate networks (20-20) on yearly TerraClimate
  `tavg_c`, `prec_mm`, `cwd_z`, recruitment limited by conspecific adult density,
  growth and mortality scored per tree. Data version `conus-annual-v1`.

Holdout metrics: see the model card in the FINN-zoo and `calibration/<variant>/report.md`
(filled in when the production fits have run).

* Repository created from the FINN model template.
