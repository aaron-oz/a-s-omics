#!/bin/bash
## Control B (R2.8) coordinate-shuffle null on klone's ckpt (idle-cycle) partition.
## One array task = one permutation. Task 0 = REAL reference; tasks 1..N = shuffles.
## Good-citizen: ckpt (preemptible idle cores), --requeue, throttled concurrency,
## summary-only output in home, transients on node-local $TMPDIR.
#SBATCH --job-name=ctrlB
#SBATCH --account=ckpt-csde          # CONFIRM with `hyakalloc`; ckpt uses your ckpt association
#SBATCH --partition=ckpt             # try ckpt / ckpt-all if ckpt is full
#SBATCH --array=0-49%20              # perm 0 (real) + 49 shuffles; max 20 running at once
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16           # mclapply parallel fits per task
#SBATCH --mem=64G                    # 16 forked INLA workers ~ 30-40G; headroom for safety
#SBATCH --time=8:00:00               # each task ~2.2h at 16 cores; buffer for a requeue after preemption
#SBATCH --requeue                    # ckpt: cleanly restart if preempted (task is idempotent)
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

set -euo pipefail
module load apptainer/1.1.5
export ASOMICS_ROOT="$HOME/a-s-omics"
export OUT="$ASOMICS_ROOT/data-outputs/mouse-embryo-raw/controls/control-B-summary.csv"
export NSAMP=1000                          # posterior draws for the median field
SIF="$HOME/containers/geospatial.sif"      # R 4.4.2 + geo stack (built once)
RLIB="$HOME/Rlib-ctrlB"                     # INLA/inlabru/fmesher/RcppHNSW/elsa (installed once)

echo "task $SLURM_ARRAY_TASK_ID on $(hostname); TMPDIR=${TMPDIR:-unset}; cpus=$SLURM_CPUS_PER_TASK"
# apptainer inherits SLURM_* env; pass the R library + our vars explicitly.
apptainer exec \
  --env R_LIBS_USER="$RLIB" \
  --env ASOMICS_ROOT="$ASOMICS_ROOT" --env OUT="$OUT" --env NSAMP="$NSAMP" \
  "$SIF" \
  Rscript "$ASOMICS_ROOT/code-spatial-smoothing/controls/control-B-klone.R"
