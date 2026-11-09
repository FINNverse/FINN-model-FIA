#!/bin/bash
# ===========================================================================
# FINN-model-FIA | run ONE calibration script on a SLURM cluster
#   sbatch calibration/hpc/run_fit.sbatch calibration/<variant>/02_fit.R
#   sbatch --dependency=afterok:<jobid> calibration/hpc/run_fit.sbatch calibration/<variant>/03_validation.R
#
# The wrapper only allocates resources and calls Rscript from the repository
# root. Every scientific choice (data, model, epochs, losses, threads) is in
# the R script itself. The scripts are single-threaded by design (they share
# nodes), so one core is enough; memory covers the 1500-site simulation of
# the validation. Expected location of the data: ../../FINN-zoo-data/fia
# relative to the repository root (see calibration/README.md).
#
# R libraries: the scripts put ~/R/lib-finn-0848b60 (FINN 0.2.0 at commit
# 0848b60) first on the search path themselves; torch and the other packages
# come from the libraries named in R_LIBS_USER below. Adapt the module and
# library lines to the cluster at hand.
# ===========================================================================
#SBATCH --job-name=finn-fia
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
SCRIPT="${1:?usage: sbatch calibration/hpc/run_fit.sbatch calibration/<variant>/<script>.R}"
module purge 2>/dev/null || true
module load math/R/4.4.2 2>/dev/null || true
export R_LIBS_USER="$HOME/R/site-library-torch017:$HOME/R/site-library"
export LD_LIBRARY_PATH="$HOME/R/site-library-torch017/torch/lib:${LD_LIBRARY_PATH:-}"
export TORCH_INSTALL=0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "[job] $SCRIPT on $(hostname), $(date)"
Rscript "$SCRIPT"
echo "[job] done $(date)"
