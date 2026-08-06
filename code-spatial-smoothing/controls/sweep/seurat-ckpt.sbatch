#!/bin/bash
## Seurat replication of the random-order sweep, on klone ckpt.
## Task 0 = ordered arm, task 1 = reverse, tasks 2-20 = 19 random draws.
## Three depths per task. Sup7's own timings were 1-40 min per evaluation, so the
## time limit is generous; each task rewrites its own output file.
#SBATCH --job-name=asomics-seurat
#SBATCH --account=ckpt-csde
#SBATCH --partition=ckpt
#SBATCH --array=0-20%21
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --output=logs/seurat-%A_%a.out
#SBATCH --error=logs/seurat-%A_%a.err

set -eo pipefail   # not -u: lmod init references unset LD_LIBRARY_PATH
module load apptainer/1.1.5

export ASOMICS_SWEEP_ROOT="${ASOMICS_SWEEP_ROOT:-$HOME/a-s-omics/sweep}"
cd "$ASOMICS_SWEEP_ROOT"
mkdir -p logs results-seurat

SIF="$HOME/containers/geospatial.sif"
RLIB="$HOME/Rlib-ctrlB"     # Seurat 5.5.1 + RcppAnnoy live here, from Control B

echo "task ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
apptainer exec \
  --env R_LIBS_USER="$RLIB" \
  --env ASOMICS_SWEEP_ROOT="$ASOMICS_SWEEP_ROOT" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  "$SIF" \
  Rscript "$ASOMICS_SWEEP_ROOT/sweep-seurat.R"
echo "task ${SLURM_ARRAY_TASK_ID} finished at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
