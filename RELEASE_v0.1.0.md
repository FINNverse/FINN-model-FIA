# Release checklist: FINNmodelFIA v0.1.0

Prepared 2026-11-10. Nothing here has been executed; tag, push and release only
after the go. All commands from the repository root.

## Commit to tag

`216c5fb` on `main` ("WorldClim licence stated correctly ..."), which follows the
release commit `5733369` ("FINNmodelFIA 0.1.0: version, FINN pinned at 0848b60,
...") and corrects the WorldClim licence wording in README, NEWS, PROVENANCE and
the release notes (WorldClim is not CC BY-SA 4.0: free for academic and other
non-commercial use, rasters not redistributed). It is the v0.1 package (`R/`, `inst/model/`, `DESCRIPTION` 0.1.0, vignette);
`calibration/` (incl. the v0.2 work in `calibration/v02/`), `hpc/` and this
file are excluded from the build by `.Rbuildignore`. `R/definition.R` also
carries the v0.2 variant definitions `oregon_v02` / `conus_v02`; they are not in
`weights.json`, so `load_model()` cannot select them, and they do not affect
the v0.1 variants (checked: both v0.1 parameter files load and simulate).

## Assets (outside git), sha256 verified 2026-11-10 with `shasum -a 256`

| asset | local path (relative to this repo) | sha256 | size |
|---|---|---|---|
| parameters_oregon.pt | ../../FINN-zoo-data/fia/fits/oregon/parameters_oregon.pt | 0d29d459a7ae10e692550c1ce64d9a94d42d8f03b7a971108ae12c57af32d815 | 3,971 B |
| model_oregon.pt | ../../FINN-zoo-data/fia/fits/oregon/model_oregon.pt | 503b7ecd37ece8cb4144f6db5dad3e9d5c5b0e18af9cf2ef9cd3783b8320092c | 4.7 MB |
| parameters_conus.pt | ../../FINN-zoo-data/fia/fits/conus/parameters_conus.pt | d40ce38f94bc307fd17fc172cf12d8c29f5b95415eb08e711d1f1a2838f2a9a9 | 17,042 B |
| model_conus.pt | ../../FINN-zoo-data/fia/fits/conus/model_conus.pt | 7c7d41915e64563cf0e2545c94c31c8cf7c0408fc0c27adee84b520f5ae38665 | 25.7 MB |

The same values are in `inst/model/weights.json` (`sha256`, `model_sha256`);
`load_model()` refuses a file whose sha256 differs.

## Checks done (2026-11-10, one thread, FINN 0.2.0 = 0848b60, torch 0.15.1, R 4.5.0)

- `R CMD build` with the vignette (chunks not evaluated at build time) and
  `R CMD INSTALL` into a temporary library: OK.
- `load_model("oregon", path = <parameters_oregon.pt>)` and
  `load_model("conus", path = <parameters_conus.pt>)`: sha256 match, model
  rebuilt (9 and 26 parameter tensors), 3-year bare-ground simulation with the
  pinned batchsize runs. `model_<variant>.pt` sha256 match.
- Vignette renders.

## Steps

```sh
# 0. sanity: the sha256 of the files to upload equal weights.json
shasum -a 256 ../../FINN-zoo-data/fia/fits/oregon/parameters_oregon.pt ../../FINN-zoo-data/fia/fits/oregon/model_oregon.pt \
               ../../FINN-zoo-data/fia/fits/conus/parameters_conus.pt  ../../FINN-zoo-data/fia/fits/conus/model_conus.pt

# 1. tag the release commit
git tag -a v0.1.0 216c5fb -m "FINNmodelFIA 0.1.0: first fitted release (variants oregon and conus), FINN 0.2.0 commit 0848b60"

# 2. push main and the tag
git push origin main
git push origin v0.1.0

# 3. create the GitHub release with the four assets
gh release create v0.1.0 --repo FINNverse/FINN-model-FIA \
  --title "FINNmodelFIA 0.1.0" \
  --notes-file RELEASE_NOTES_v0.1.0.md \
  ../../FINN-zoo-data/fia/fits/oregon/parameters_oregon.pt \
  ../../FINN-zoo-data/fia/fits/oregon/model_oregon.pt \
  ../../FINN-zoo-data/fia/fits/conus/parameters_conus.pt \
  ../../FINN-zoo-data/fia/fits/conus/model_conus.pt

# 4. verify the uploaded assets byte for byte
tmp=$(mktemp -d) && gh release download v0.1.0 --repo FINNverse/FINN-model-FIA -D "$tmp" && shasum -a 256 "$tmp"/*.pt

# 5. (user's go) make the repository public
gh repo edit FINNverse/FINN-model-FIA --visibility public --accept-visibility-change-consequences
```

## After the release

```r
remotes::install_github("FINNverse/FINN-model-FIA@v0.1.0")
library(FINNmodelFIA); torch::torch_set_num_threads(1L)
m <- load_model("conus")      # downloads from the release URL in weights.json and checks the sha256
```

Then the FINN-zoo card (`FIA/model-card.qmd`) cites tag `v0.1.0` and these
sha256 values.

## Release notes

The text of `RELEASE_NOTES_v0.1.0.md` (used by `--notes-file` above).
