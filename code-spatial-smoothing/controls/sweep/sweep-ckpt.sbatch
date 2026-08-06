#!/bin/bash
#SBATCH --job-name=asomics-sweep
#SBATCH --account=ckpt-csde
#SBATCH --partition=ckpt
#SBATCH --requeue
#SBATCH --array=0-101%40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/sweep-%A_%a.out
#SBATCH --error=logs/sweep-%A_%a.err

# Random-order feature-depth sweep (see SCOPE.md).
#
# 102 tasks:
#   0        ordered arm over the seven design depths
#   1        reverse arm over the same
#   2        ordered arm over all ~30 published depths (reference curve)
#   3-101    99 random draws, seven depths each
#
# Each task is single-threaded and takes about 2.5 minutes, except task 2 which is
# about 13 minutes. Every task rewrites its own output file from scratch, so preemption
# and requeue are safe with no partial-state handling.

set -euo pipefail

export ASOMICS_SWEEP_ROOT="${ASOMICS_SWEEP_ROOT:-$HOME/a-s-omics/sweep}"
cd "$ASOMICS_SWEEP_ROOT"
mkdir -p logs results

# keep BLAS single-threaded: the clustering is deterministic only at one thread, and the
# array already provides the parallelism
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# ---- R environment -------------------------------------------------------------------
# Needs only Matrix, RcppHNSW, igraph, elsa, raster. No INLA, no Seurat, unlike Control B.
# Set the right line for the current klone build before submitting.
module load r/4.4.0 2>/dev/null || module load R 2>/dev/null || true

echo "task ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
Rscript sweep-run.R
echo "task ${SLURM_ARRAY_TASK_ID} finished at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
