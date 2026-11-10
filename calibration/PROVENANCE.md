# Provenance of the calibration data

Every raw layer the calibration chain reads, with source URL, licence, date pulled and
sha256. Paths are relative to the data folder `../../FINN-zoo-data/fia` (outside git,
two levels above the repository root; see README.md in this folder). Both variants
start from the public FIADB dump; the whole chain from that dump to the prepared sets
is in this folder (`conus/01..09`, `oregon/01`) and has been checked to reproduce the
tables of the FINN-fia analysis pipeline byte for byte (`expected_preparation.csv`).

## 1. FIA Database, entire US, CSV dump (both variants)

Source: USDA Forest Service, Forest Inventory and Analysis DataMart,
<https://apps.fs.usda.gov/fia/datamart/datamart.html> ("Download all states combined",
CSV). Licence: US Government work, public domain (CC0 1.0 as stated on the DataMart).
Pulled 2025-06-23/24. Expected under `raw/fiadb/` (or the FINN-fia data folder).

| file | size (bytes) | sha256 |
|---|---|---|
| `CSV_FIADB_ENTIRE/ENTIRE_SURVEY.csv` | 311,159 | `1070dd193b0f446c500515b9e4a4f4176f95cfe03d04f16c3c55bb468e398097` |
| `CSV_FIADB_ENTIRE/ENTIRE_PLOT.csv` | 451,669,663 | `1dd10b3f98bbdaf0a984e058d20065de0b94dd5f67d694126b6797a4b141776f` |
| `CSV_FIADB_ENTIRE/ENTIRE_SUBPLOT.csv` | 1,063,678,146 | `3035e03d13581746c8086a20a57bdf20a039989b204cdec033207ca32361ed18` |
| `CSV_FIADB_ENTIRE/ENTIRE_COND.csv` | 754,235,850 | `fb17c992c71d90ba1e5e1fd98748eae08c4d09571e202dfb53a3a7f277ef65d5` |
| `CSV_FIADB_ENTIRE/ENTIRE_TREE.csv` | 12,899,145,944 | `7427bb534925713bf5696a08805dfbb733dd53d0ba95051de079ef98cb8a4d01` |
| `FIADB_REFERENCE/REF_SPECIES.csv` | 592,260 | `179124c82da156b63ec766a92f6fded3e5e9d24da060f974588628212e2dee2f` |

Read by `conus/01_fiadb_metadata.R` (all six), `conus/03_extract_trees.R` (TREE,
REF_SPECIES) and `oregon/01_prepare.R` (PLOT, for the Oregon plot coordinates).

## 2. WorldClim 2.1, 1970-2000 normals, 30 arc-seconds (both variants)

Source: <https://www.worldclim.org/data/worldclim21.html> (bioclimatic variables,
30 s; the bundle `wc2.1_30s_bio.zip`, 10,401,038,857 bytes). Licence: WorldClim's own
terms, free for academic and other non-commercial use, no redistribution without
permission (https://www.worldclim.org/about.html); the rasters are not redistributed,
only values extracted at the plot locations are included. Attribution Fick, S.E. &
Hijmans, R.J. (2017) *Int. J. Climatol.* 37, 4302-4315, doi:10.1002/joc.5086. Rasters
dated 2019-02-28; bundle pulled 2025-09-09. Expected under `raw/worldclim/` (or a
`worldclim` folder next to the models folder).

| file | size (bytes) | sha256 |
|---|---|---|
| `wc2.1_30s_bio/wc2.1_30s_bio_1.tif` (BIO1, annual mean temperature) | 753,795,566 | `382c991cea4a20b42e345fcd5a897788feee05bc0d613dec032ff0cad1f4c375` |
| `wc2.1_30s_bio/wc2.1_30s_bio_12.tif` (BIO12, annual precipitation) | 487,101,349 | `6df0ef2f57f900733a569fcc4c8debd27b490978c5d82c2970485f0e66abbf28` |

Read by `conus/04_attach_worldclim.R`. The Oregon export (section 4) carries the same
two variables, plus BIO5, BIO6, BIO15, BIO18, extracted upstream from these rasters.

## 3. TerraClimate, monthly, 4 km, 1958-2023 (variant conus)

Source: Abatzoglou, J.T. et al. (2018) *Sci. Data* 5, 170191; served by the Northwest
Knowledge Network OPeNDAP service
`https://thredds.northwestknowledge.net/thredds/dodsC/TERRACLIMATE_ALL/data/TerraClimate_<var>_<year>.nc`
(variables `def`, `vpd` 1958-2023; `tmax`, `tmin`, `ppt` 2000-2024). Licence: public
domain (CC0). This is a live API: the plot-level monthly series pulled on 2026-09-08 for
the 33,320 cohort plots are the reproducible layer; `conus/05` and `conus/06` reuse them
when they cover the cohort exactly and pull again otherwise.

| file | size (bytes) | sha256 |
|---|---|---|
| `intermediate/conus/drought_monthly_cache.rds` (def, vpd; 33,320 plots x 12 months x 1958-2023) | 64,232,816 | `8b29378cb6562518818459551b28fe36badc008987e7c78baaadbc09cf7ff276` |
| `intermediate/conus/annual_climate_cache.rds` (tmax, tmin, ppt; 33,320 plots x 12 months x 2000-2024) | 49,264,966 | `5082d2fdb4dba89dc54391933909e0c54c6344a7f7617d462e801d7b56384edf` |

## 4. Oregon export (variant oregon)

The Oregon tree list and climate that ship with the FINN package sources
(`FINN/data-raw/`, commit 15922a1, 2026-06-24), written by the FINN-fia pipeline
(`scripts/03_attach_environment.R`, `07_prepare_finn_inputs.R`, interval mode) from the
FIADB dump of section 1 and the WorldClim rasters of section 2. Derived from
public-domain FIA data and WorldClim values extracted at the plot locations (free for
academic and other non-commercial use); distributed with FINN under GPL-3.
Copied 2026-11-09 to `raw/oregon-finn-data-raw/`.

| file | sha256 |
|---|---|
| `full_tree_dt.csv` (775 plots, 4 subplots, inventories 2001-03 / 2011-13 / 2021-24) | `4a755925711c669dce7c1c242b6d96e08895afe0206abe284ccd44b0a02c9fdd` |
| `env_unscaled_dt.csv` (WorldClim normals per plot: temp, tempmax, tempmin, prec, precseas, precwarmq) | `451b463e47356d1ac05c4818e4da9c187df4419d2a703fca05c04a281696c861` |
| `species_dt.csv` (27 species + other) | `3449ea1c24a41d5bce109baeb18774616e3d27a60c4c143cb53af98e9f0fa958` |
| `oregon_plot_coords.csv` (Oregon rows of ENTIRE_PLOT.csv: COUNTYCD, PLOT, LAT, LON; written by `oregon/01_prepare.R`, 2026-11-09) | `43b9d2180ce32fe82747472cbefb2edfe532234d5dee6808eb57fe3b0be3f6c9` |

## Reproduction check (2026-11-09)

`conus/01..09` were run on the laptop against the layers above (TerraClimate from the
caches). The resulting `prepared/annual-v5/{train1500,holdout_blocked}` are byte-identical
to the FINN-fia tables of 2026-09-10 (obs_dt, tree_obs_dt, env_dt, env_unscaled_dt,
init_cohort_dt, species_dt, siteID_dt, site_calendar, site_blocks), and
`expected_preparation.csv` is unchanged. One constant mattered: the acre-to-hectare
factor of the FINN-fia pipeline (0.404686) is kept in `conus/08`; with the exact factor
0.40468564224 the recruitment column differs in the 7th significant digit.

## Redistribution

FIA plot data are public domain; the FIADB coordinates are already perturbed by the
Forest Service, and the FINN tables carry no personal data. The example sites in
`inst/extdata` are six holdout plots per variant (environment and initial tree list)
and may be redistributed. Climate normals: WorldClim 2.1 (Fick & Hijmans 2017,
doi:10.1002/joc.5086), free for academic and other non-commercial use; the rasters are
not redistributed, only values extracted at the plot locations are included.
