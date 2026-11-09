# Provenance of the calibration data

Every raw layer the calibration chain reads, with source, licence, date pulled and
sha256. Paths are relative to the data folder `../../FINN-zoo-data/fia` (outside git,
two levels above the repository root; see README.md in this folder).

## Common root: the FIADB CSV dump

| layer | source | licence | pulled | file (sha256) |
|---|---|---|---|---|
| FIA Database, entire US, CSV | USDA Forest Service, FIA DataMart, https://apps.fs.usda.gov/fia/datamart/datamart.html ("Download all states combined", CSV) | US Government work, public domain (CC0 1.0 as stated on the DataMart) | 2025-06-23/24 | `ENTIRE_PLOT.csv` 451,669,663 B, sha256 `1dd10b3f98bbdaf0a984e058d20065de0b94dd5f67d694126b6797a4b141776f`; `ENTIRE_TREE.csv` 12,899,145,944 B; `ENTIRE_COND.csv`, `ENTIRE_SUBPLOT.csv`, `ENTIRE_SUBP_COND.csv`, `ENTIRE_SURVEY.csv`, `REF_SPECIES.csv` (the FINN-fia pipeline reads these; only `ENTIRE_PLOT.csv` is read directly here) |

The dump lives in the FINN-fia analysis repository (`data/CSV_FIADB_ENTIRE/`); the
coordinates step of the Oregon variant reads `ENTIRE_PLOT.csv` from there
(`PARAMS$fiadb_plot_csv` in `oregon/01_prepare.R`) and caches the Oregon rows.

## Variant `oregon`

| layer | source | licence | pulled | file (sha256) |
|---|---|---|---|---|
| Oregon tree list, three inventories | `FINN/data-raw/full_tree_dt.csv` (FINN package sources, commit 15922a1, 2026-06-24), exported by FINN-fia `scripts/03_attach_environment.R` + `07_prepare_finn_inputs.R` from the FIADB dump | derived from public-domain FIA data; GPL-3 with FINN | copied 2026-11-09 | `raw/oregon-finn-data-raw/full_tree_dt.csv` sha256 `4a755925711c669dce7c1c242b6d96e08895afe0206abe284ccd44b0a02c9fdd` |
| plot climate (WorldClim 2.1 normals 1970-2000, 30 arc-sec: bio1, bio5, bio6, bio12, bio15, bio18) | same export, `env_unscaled_dt.csv`; WorldClim 2.1, https://www.worldclim.org/data/worldclim21.html | WorldClim: CC BY-SA 4.0 (attribution: Fick & Hijmans 2017) | copied 2026-11-09 | `raw/oregon-finn-data-raw/env_unscaled_dt.csv` sha256 `451b463e47356d1ac05c4818e4da9c187df4419d2a703fca05c04a281696c861` |
| species codes | same export, `species_dt.csv` | as above | copied 2026-11-09 | `raw/oregon-finn-data-raw/species_dt.csv` sha256 `3449ea1c24a41d5bce109baeb18774616e3d27a60c4c143cb53af98e9f0fa958` |
| plot coordinates (public, perturbed LAT/LON) | FIADB `ENTIRE_PLOT.csv`, Oregon rows, cached by `oregon/01_prepare.R` | public domain | 2026-11-09 | `raw/oregon-finn-data-raw/oregon_plot_coords.csv` sha256 `43b9d2180ce32fe82747472cbefb2edfe532234d5dee6808eb57fe3b0be3f6c9` |

## Variant `conus`

The conus variant starts from the prepared set `annual-v5` of the FINN-fia analysis
pipeline (private repository FINNverse/FINN-fia, `outputs/prepared/annual-v5`, written
2026-09-10 by `scripts/07_prepare_finn_inputs.R` at FINN-fia commit 58fcf3e). Its
preparation from the public sources (FIADB dump above; TerraClimate
https://www.climatologylab.org/terraclimate.html, CC0/public domain, monthly 1958-2023,
summarised to yearly plot values by FINN-fia `scripts/03c_attach_drought.R` and
`03d_attach_annual_climate.R`; WorldClim 2.1 as above) is documented in the appended
FINN-fia PREPARATION.md (`prepared/conus-annual-v1/PREPARATION.md`) and is NOT
repeated by this repository. Reproducing the conus tables from scratch therefore
needs the FINN-fia code; this is an open point for the first public release.

| layer | source | licence | pulled | file (sha256) |
|---|---|---|---|---|
| train1500 observations | FINN-fia annual-v5 | derived from public-domain FIA + TerraClimate (CC0) | copied 2026-11-09 | `prepared/annual-v5/train1500/obs_dt.csv` `c3e34e733911eff0b1b2c6adecb7131f920351ac0212ce4152f71c366f73d017` |
| train1500 tree intervals | FINN-fia annual-v5 | as above | copied 2026-11-09 | `prepared/annual-v5/train1500/tree_obs_dt.csv` `877350d1e41d46f2c7f93b7397f46c33a4e7107c4f7f3275f7f6248bde8f9417` |
| train1500 environment (natural units) | FINN-fia annual-v5 | as above | copied 2026-11-09 | `prepared/annual-v5/train1500/env_unscaled_dt.csv` `bef6e21d8f21e3edc2e68f4592bba2059f8b9594216a04acd6b5c83fe26245c9` |
| train1500 initial trees | FINN-fia annual-v5 | as above | copied 2026-11-09 | `prepared/annual-v5/train1500/init_cohort_dt.csv` `466a9b40d1808a34d3e9ef3139dae284be704cc4fc7579b6c42f74b0da3d4320` |
| holdout_blocked observations | FINN-fia annual-v5 | as above | copied 2026-11-09 | `prepared/annual-v5/holdout_blocked/obs_dt.csv` `6b50227e19a45e82291d75e101530b42db8d8004ad3ae96c2e893d783d6b3bab` |
| holdout_blocked tree intervals | FINN-fia annual-v5 | as above | copied 2026-11-09 | `prepared/annual-v5/holdout_blocked/tree_obs_dt.csv` `016859913f75a9d924eac59d0a886b746ec061c7770cb0494c9487b2cab1d718` |
| holdout_blocked environment | FINN-fia annual-v5 | as above | copied 2026-11-09 | `prepared/annual-v5/holdout_blocked/env_unscaled_dt.csv` `c78a58544a774955c0797643dc4aa036ee09463fd5b621824416a091848aa134` |
| holdout_blocked initial trees | FINN-fia annual-v5 | as above | copied 2026-11-09 | `prepared/annual-v5/holdout_blocked/init_cohort_dt.csv` `4ec6c72e78bb12d9828d76f72cb673a7bce89dbf5f78b6903eb288e7caadc944` |

## Redistribution

FIA plot data are public domain; the plot coordinates in FIADB are already perturbed
by the Forest Service, and the FINN tables carry no personal data. The example sites
in `inst/extdata` are six holdout plots per variant (environment and initial tree
list) and may be redistributed. WorldClim-derived values require attribution
(Fick & Hijmans 2017, CC BY-SA 4.0).
