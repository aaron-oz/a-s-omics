#!/bin/bash
## Random-order feature-depth sweep on klone's ckpt (idle-cycle) partition. See SCOPE.md.
##
## 102 array tasks:
##   0        ordered arm over the seven design depths
##   1        reverse arm over the same
##   2        ordered arm over all ~30 published depths (the reference curve)
##   3-101    99 random draws, seven depths each
##
## Each task is single-threaded, about 2.5 min, except task 2 at about 13 min. Every task
## rewrites its own output CSV from scratch, so --requeue after preemption needs no
## partial-state handling.
##
## Reuses the container and R library built for Control B: no new environment setup.
#SBATCH --job-name=asomics-sweep
#SBATCH --account=ckpt-csde          # CONFIRM with `hyakalloc`
#SBATCH --partition=ckpt             # ckpt-all widens the node pool if ckpt is busy
#SBATCH --array=0-101%20  # throttled: ckpt is shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G          # measured MaxRSS 3.7G; 2x headroom
#SBATCH --time=01:00:00              # ~2.5 min typical; buffer for a requeue after preemption
#SBATCH --requeue
#SBATCH --output=logs/sweep-%A_%a.out
#SBATCH --error=logs/sweep-%A_%a.err

set -euo pipefail
module load apptainer/1.1.5

export ASOMICS_SWEEP_ROOT="${ASOMICS_SWEEP_ROOT:-$HOME/a-s-omics/sweep}"
cd "$ASOMICS_SWEEP_ROOT"
mkdir -p logs results

SIF="$HOME/containers/geospatial.sif"   # R 4.4.2 + geo stack, built for Control B
RLIB="$HOME/Rlib-ctrlB"                 # RcppHNSW, elsa; igraph/raster/Matrix are in the image

# The clustering is deterministic only single-threaded, and the array supplies the
# parallelism, so pin every BLAS/OpenMP layer to one thread.
echo "task ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
apptainer exec \
  --env R_LIBS_USER="$RLIB" \
  --env ASOMICS_SWEEP_ROOT="$ASOMICS_SWEEP_ROOT" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  "$SIF" \
  Rscript "$ASOMICS_SWEEP_ROOT/sweep-run.R"
echo "task ${SLURM_ARRAY_TASK_ID} finished at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
