#!/bin/bash
# ===========================================================================
# FINN-model-FIA | run ONE v0.2 pilot (or any calibration script) on a SLURM cluster
#   sbatch hpc/run_pilot.sbatch calibration/v02/pilot_oregon_O0_v01.R
#
# Resource-only wrapper: one CPU, one thread, 2 h, 8 GB; calls Rscript from the
# repository root. Every scientific choice is inside the R script. The scripts put
# ~/R/lib-finn-0848b60 first on .libPaths() themselves and fall back to
# torch::optim_adam where optim_ignite_adam is missing (torch 0.13). Data are
# expected under ../../FINN-zoo-data/fia relative to the repository root
# (calibration/v02/CLUSTER.md lists what to sync).
# ===========================================================================
#SBATCH --job-name=fia-pilot
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
mkdir -p logs
SCRIPT="${1:?usage: sbatch hpc/run_pilot.sbatch calibration/v02/<script>.R}"
module purge 2>/dev/null || true
module load math/R/4.4.2 2>/dev/null || true
export R_LIBS_USER="$HOME/R/site-library-torch017:$HOME/R/site-library"
export LD_LIBRARY_PATH="$HOME/R/site-library-torch017/torch/lib:${LD_LIBRARY_PATH:-}"
export TORCH_INSTALL=0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "[job] $SCRIPT on $(hostname), $(date)"
Rscript "$SCRIPT"
echo "[job] done $(date)"
