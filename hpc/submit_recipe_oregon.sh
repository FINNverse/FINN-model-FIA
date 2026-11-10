#!/bin/bash
# ===========================================================================
# FINN-model-FIA | submit the zoo v0.2 recipe chain of the OREGON variant on the cluster
#   bash hpc/submit_recipe_oregon.sh            (from the repository root)
#
# Three jobs, each one self-contained R script run by hpc/run_pilot.sbatch
# (resources only; every scientific choice is in the script):
#   stage A  calibration/v02/recipe/oregon/fit_stageA.R   demography, observed recruitment
#   stage B  calibration/v02/recipe/oregon/fit_stageB.R   regeneration (afterok: stage A)
#   report   calibration/v02/recipe/oregon/validation.R   holdout report (afterok: stage B)
# Time limits from the measured cost (laptop M-series, one thread, 2026-11-10):
# stage A 0.4 s/epoch on the first 200 training plots -> ~1 s/epoch on the 450
# plots, 300 epochs ~5 min; stage B the same per epoch, 150 epochs ~3 min. The
# cluster's cores ran the v0.2 pilots ~3x slower than the laptop (Oregon pilots
# 9-12 min for 300 epochs), hence 2 h and 1 h below. Validation: 561 plots over two
# steps plus six 50-step bare-ground runs, 4 threads.
# ===========================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
A=$(sbatch --parsable --job-name=zoo-fia-oregon-A --time=02:00:00 --mem=8G hpc/run_pilot.sbatch calibration/v02/recipe/oregon/fit_stageA.R)
B=$(sbatch --parsable --job-name=zoo-fia-oregon-B --time=01:00:00 --mem=8G --dependency=afterok:"$A" hpc/run_pilot.sbatch calibration/v02/recipe/oregon/fit_stageB.R)
V=$(sbatch --parsable --job-name=zoo-fia-oregon-V --time=01:00:00 --mem=8G --cpus-per-task=4 --dependency=afterok:"$B" hpc/run_pilot.sbatch calibration/v02/recipe/oregon/validation.R)
echo "submitted: stage A $A -> stage B $B (afterok) -> validation $V (afterok); logs in logs/fia-pilot_<jobid>.out|err"
