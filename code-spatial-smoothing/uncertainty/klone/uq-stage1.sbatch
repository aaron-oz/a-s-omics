#!/bin/bash
## Stage 1: refit genes and write int16 joint draws. Preemptible ckpt with requeue; tasks are
## idempotent (completed genes are skipped). One CPU: INLA runs single-threaded.
## SIZING: --mem and --time below are for the ONE measurement task; the array is resubmitted with
## values sized from sacct (MaxRSS x 2, Elapsed + buffer) and the measured numbers recorded here.
#SBATCH --job-name=uq-s1
#SBATCH --account=ckpt-csde
#SBATCH --partition=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:45:00
#SBATCH --output=logs/uq-s1-%A_%a.out
#SBATCH --error=logs/uq-s1-%A_%a.err
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export ROOT=$HOME/a-s-omics/uq DRAWS=${DRAWS:-/gscratch/csde/azimmer/uq-draws} OUT=$HOME/a-s-omics/uq/out BATCH=${BATCH:-10}
mkdir -p "$DRAWS" "$OUT"
module load apptainer/1.1.5
apptainer exec --bind /gscratch --env R_LIBS_USER=$HOME/Rlib-ctrlB $HOME/containers/geospatial.sif \
  Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); source(file.path(Sys.getenv("ROOT"), "uq-stage1.R"))'
