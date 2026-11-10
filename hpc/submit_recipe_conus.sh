#!/bin/bash
# ===========================================================================
# FINN-model-FIA | submit the zoo v0.2 recipe chain of the CONUS variant on the cluster
#   bash hpc/submit_recipe_conus.sh            (from the repository root)
#
# Three jobs, each one self-contained R script run by hpc/run_pilot.sbatch
# (resources only; every scientific choice is in the script):
#   stage A  calibration/v02/recipe/conus/fit_stageA.R   demography, observed recruitment
#   stage B  calibration/v02/recipe/conus/fit_stageB.R   regeneration (afterok: stage A)
#   report   calibration/v02/recipe/conus/validation.R   holdout report (afterok: stage B)
# Time limits from the measured cost (laptop M-series, one thread, 2026-11-10):
# stage A 2.7 s/epoch on the first 200 training plots -> ~20 s/epoch on the 1,500
# plots, 300 epochs ~1.7 h; stage B 2.4 s/epoch on 200 -> ~18 s/epoch, 150 epochs
# ~0.8 h. The cluster's cores ran the v0.2 pilots ~3x slower than the laptop, and
# each fit script ends with a full simulation of the training set, hence 12 h and
# 6 h below (a 2x margin on top of the factor 3). Validation: 2,500 plots over 23
# years plus six 500-year bare-ground runs, 4 threads.
# ===========================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
A=$(sbatch --parsable --job-name=zoo-fia-conus-A --time=12:00:00 --mem=16G hpc/run_pilot.sbatch calibration/v02/recipe/conus/fit_stageA.R)
B=$(sbatch --parsable --job-name=zoo-fia-conus-B --time=06:00:00 --mem=16G --dependency=afterok:"$A" hpc/run_pilot.sbatch calibration/v02/recipe/conus/fit_stageB.R)
V=$(sbatch --parsable --job-name=zoo-fia-conus-V --time=03:00:00 --mem=16G --cpus-per-task=4 --dependency=afterok:"$B" hpc/run_pilot.sbatch calibration/v02/recipe/conus/validation.R)
echo "submitted: stage A $A -> stage B $B (afterok) -> validation $V (afterok); logs in logs/fia-pilot_<jobid>.out|err"
