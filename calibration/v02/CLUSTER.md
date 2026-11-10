# v0.2 pilots on the cluster: what to sync and how to run

The eight pilot scripts are self-contained (`Rscript <script>` from the repository root,
one thread, each <= 45 min on the laptop; a 2 h limit covers the cluster's slower
`optim_adam` path). They read data from `../../FINN-zoo-data/fia` relative to the
repository root, so the cluster needs the same two-level layout:

```
<base>/FINN-zoo-models/FINN-model-FIA/      this repository (branch state of the commit named in the notes)
<base>/FINN-zoo-data/fia/                   data, outside git
```

All pilots run on the v0.2 DATA (agreed rare-species rule, 2026-11-10): Oregon
`oregon-10yr-v2` (S = 18 species, < 2 % out-of-S basal area at the interval-start
censuses, 450 train / 111 holdout plots), CONUS `conus-annual-v2` (S = 36 species,
same rule, complete ingrowth target; the pilots draw 300 train / 200 holdout plots with
fixed seeds). O0 / C0 are the v0.1 settings on these data.

## Files of the repository that the pilots use
- `calibration/v02/pilot_eval.R` (shared helpers: evaluation, stratification by the
  end-census out-of-S share, theta-by-name mapping, recruitment cap; sourced)
- `R/definition.R` (sourced by O2, O3, C1, C2, C3 for the v0.1 theta via `build_model()`)
- `calibration/oregon/env_scaling.csv`, `calibration/conus/env_scaling.csv` (v0.1 scaling, needed by `build_model()`)
- the pilot scripts: `calibration/v02/pilot_oregon_O0_v01.R`, `pilot_oregon_O1_cap_cosine.R`,
  `pilot_oregon_O2_v2_cap_theta.R`, `pilot_oregon_O3_v2_staged.R`, `pilot_conus_C0_v01.R`,
  `pilot_conus_C1_v2_cap_theta_cosine.R`, `pilot_conus_C2_v2_staged.R`
- `calibration/v02/pilot_conus_C3_floor.R` (the mortality-floor arm: needs the branch build
  `zoo-v02-mortfloor` installed as `~/R/lib-finn-zoo-v02` on the cluster; the script puts
  that library first and stops if `createProcess(mortality_floor=)` is missing)
- `hpc/run_pilot.sbatch` (resource-only wrapper: 1 CPU, `OMP_NUM_THREADS=1`, 8 GB, 2 h)
- `calibration/v02/compare_pilots.R` (afterwards, on the copied-back results)

## Data folders to sync (from `FINN-zoo-data/fia/`, nothing else)
| folder | used by | size |
|---|---|---|
| `prepared/oregon-10yr-v2/{train_blocked,holdout_blocked}/` | O0-O3 | 2 MB |
| `prepared/conus-annual-v2/{train,holdout}/` | C0-C3 | 19 MB |
| `prepared/oregon-10yr-v1/train_blocked/species_dt.csv`, `prepared/conus-annual-v1/train/species_dt.csv` | O2, O3, C1, C2, C3 (the v0.1 species names for the theta mapping) | 2 KB |
| `fits/oregon/parameters.pt` (4 KB), `fits/conus/parameters.pt` (17 KB) | O2, O3, C1, C2, C3 (v0.1 theta) | 21 KB |

The pilots write into `FINN-zoo-data/fia/fits/pilot_<name>.pt` (parameter lists) and
`FINN-model-FIA/calibration/v02/pilots/<name>/` (metrics.csv, summary.csv,
metrics_by_end_share.csv, predictions.csv, bare_ground.csv, loss_history.csv). Copy the
`calibration/v02/pilots/` folder back; `Rscript calibration/v02/compare_pilots.R` then
writes `PILOTS_v02.md`.

## Submission (one job per pilot, independent of each other)
```
cd <base>/FINN-zoo-models/FINN-model-FIA
for p in pilot_oregon_O0_v01 pilot_oregon_O1_cap_cosine pilot_oregon_O2_v2_cap_theta pilot_oregon_O3_v2_staged \
         pilot_conus_C0_v01 pilot_conus_C1_v2_cap_theta_cosine pilot_conus_C2_v2_staged pilot_conus_C3_floor; do
  sbatch hpc/run_pilot.sbatch calibration/v02/$p.R
done
```
Library: `~/R/lib-finn-0848b60` (FINN 0.2.0 at 0848b60; its torch 0.13 lacks
`optim_ignite_adam`, every script falls back to `torch::optim_adam`); C3 additionally
`~/R/lib-finn-zoo-v02`. Expected cost: Oregon pilots 300 epochs x ~2-4 s = 15-25 min
each (450 plots, 18 species); conus pilots 250 epochs x ~15-25 s on 300 plots with 36
species = 60-100 min each (the 36-species cohort tensors are ~1.8x the v0.1 size; if a
conus pilot exceeds 2 h, resubmit with `--time=04:00:00`); the evaluation adds 3-10 min.
Memory < 4 GB.

The laptop pilot queue is stopped (nothing runs twice).

## v0.2 production fits: the recipe chains (2026-11-10)

The production fits follow the zoo v0.2 recipe (`FINN-zoo-data/RECIPE_v02.md`), one
chain per variant, each three self-contained scripts (`Rscript <script>` from the
repository root, one thread for the fits, four for the validation):

| variant | scripts | data | fit folder (outside git) | chain |
|---|---|---|---|---|
| conus | `calibration/v02/recipe/conus/{fit_stageA.R, fit_stageB.R, validation.R}` | `prepared/conus-annual-v2/{train,holdout}` | `fits/conus-annual-v2-recipe/` | `bash hpc/submit_recipe_conus.sh` (A 12 h / B 6 h, 16 GB; V 3 h, 4 CPUs) |
| oregon | `calibration/v02/recipe/oregon/{fit_stageA.R, fit_stageB.R, validation.R}` | `prepared/oregon-10yr-v2/{train_blocked,holdout_blocked}` | `fits/oregon-10yr-v2-recipe/` | `bash hpc/submit_recipe_oregon.sh` (A 2 h / B 1 h, 8 GB; V 1 h, 4 CPUs) |

Stage A fits competition, growth and the U-shaped mortality with climate networks while
the observed recruits are fed in; stage B fits `regeneration_saturation` with the
demography and the ML dispersion frozen (all losses); validation writes `report.md` and
`validation/` into the fit folder. The chains use `--dependency=afterok`. Laptop smoke
tests (3 epochs, the first 200 training plots, one thread): conus 2.7 / 2.4 s per epoch
(A / B), Oregon 0.4 / 0.4 s; the limits above carry the cluster's factor 3 and a 2x
margin. The drafts `calibration/v02/fit_conus.R` and `fit_oregon.R` predate the recipe
and are superseded by the recipe folders.

Deviations from the reference implementation (documented in every script header):
the observed recruits are fed as a per-patch mean rate and the recruits of every stage
and simulation are drawn with FINN's `regeneration_sampler = "refined"` (a plot is
4 subplots of 0.01681 ha with one site-level env, so one observed recruit is a quarter
tree per patch, which FINN's deterministic rounding and its default "gaussian" sampler
drop); conus uses the yearly TerraClimate columns incl. `cwd_z` (the data version has no
normals) and per-interval windows for plots with several censuses; Oregon has one step
per decade and no tree table (site-level growth and mortality losses), and uses the
climate network the recipe asks for (the smoke test fits it on 200 plots). K is
trainable as the recipe says; the stage-B summary reports it against the observed
recruitment maxima so a runaway (the v0.1 Oregon fit) is visible.

Fetch the results: `rsync -av <cluster>:finn-zoo/FINN-zoo-data/fia/fits/conus-annual-v2-recipe/ ../../FINN-zoo-data/fia/fits/conus-annual-v2-recipe/`
and the same for `oregon-10yr-v2-recipe/`.
