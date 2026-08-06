#!/bin/bash
## Control A v2: exclusion-aware ligand-receptor pair-shuffle null, on klone ckpt.
## Task 0 = the true pairing; tasks 1-99 = scrambles with zero real pairs remaining.
## Each task is single-threaded, roughly 1-2 min, and rewrites its own output file.
#SBATCH --job-name=asomics-ctrlA
#SBATCH --account=ckpt-csde
#SBATCH --partition=ckpt
#SBATCH --array=0-99%40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --output=logs/ctrlA-%A_%a.out
#SBATCH --error=logs/ctrlA-%A_%a.err

set -eo pipefail   # not -u: lmod init references unset LD_LIBRARY_PATH
module load apptainer/1.1.5

export ASOMICS_CTRLA_ROOT="${ASOMICS_CTRLA_ROOT:-$HOME/a-s-omics/sweep}"
cd "$ASOMICS_CTRLA_ROOT"
mkdir -p logs results-ctrlA

SIF="$HOME/containers/geospatial.sif"
RLIB="$HOME/Rlib-ctrlB"

echo "task ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
apptainer exec \
  --env R_LIBS_USER="$RLIB" \
  --env ASOMICS_CTRLA_ROOT="$ASOMICS_CTRLA_ROOT" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  "$SIF" \
  Rscript "$ASOMICS_CTRLA_ROOT/control-A-v2.R"
echo "task ${SLURM_ARRAY_TASK_ID} finished at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
