# FINN model repositories: conventions

Each fitted FINN model the user wants to offer as a tool lives in its own
repository in the FINNverse organisation (`FINN-model-<data>`), packaged as a
small R package (`FINNmodel<Data>`). The repo is NOT tied to a manuscript: our
papers may use similar models, but the tool is versioned for users and evolves
on its own. The FINN-zoo site is the catalogue: one model card per model,
linking to the repo and its current release.

## Layout
```
DESCRIPTION, NAMESPACE, LICENSE (GPL-3), README.md, NEWS.md
R/model.R          shared across all model repos; do not change per repo
R/definition.R     build_model(variant, parameters): the model structure, written out
inst/model/weights.json   per variant: release tag, file, url, sha256, FINN commit, patches, patch_size
inst/model/species.csv    species index -> name (and any codes)
inst/model/covariates.csv env covariates: name, unit, source, how to derive for new sites
inst/extdata/      a few example sites (env + init cohorts) that may be redistributed
vignettes/quickstart.Rmd  load, simulate example sites, bare-ground run
calibration/       the self-contained pipeline that produced the release:
                   01-... download/prepare, fit.R, validation.R, PREPARATION.md, report.md
                   (house rules: r-pipeline-conventions, finn-standalone-fit)
```
`.pt` files and data are never committed (`.gitignore`).

## Weights and versions
- Saved object: `torch_save()` of a NAMED LIST of parameter tensors (not the
  whole model object), so `build_model()` can rebuild with the installed FINN.
- A release = git tag `vX.Y.Z` with `parameters.pt` (one per variant) as release
  asset; `weights.json` carries tag, url, sha256 and FINN commit.
- Semver: major = species set, covariates or model structure change (users must
  change inputs); minor = refit on more/new data with the same inputs; patch =
  code/doc fixes, same weights.
- `NEWS.md` per release: data version, FINN commit, holdout metrics.
- Repositories start private; they become public at the first release, with the
  user's go.

## Rules
- Data redistribution: only what the source license allows (example sites in
  `inst/extdata` only if allowed; BWI: 1 km INSPIRE cells only, never personal
  names). The calibration scripts download/prepare the rest from the source.
- No pushing by agents; the coordinator pushes and reports. No Claude
  attribution in commits; no host names or absolute home paths in committed files.
