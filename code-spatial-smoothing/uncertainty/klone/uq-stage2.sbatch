#!/bin/bash
## Stage 2: per pair, minimum of ligand and receptor draws, region means, slice draws. Preemptible ckpt with requeue; tasks are
## idempotent (completed genes are skipped). One CPU: INLA runs single-threaded.
## SIZING: --mem and --time below are for the ONE measurement task; the array is resubmitted with
## values sized from sacct (MaxRSS x 2, Elapsed + buffer) and the measured numbers recorded here.
#SBATCH --job-name=uq-s2
#SBATCH --account=ckpt-csde
#SBATCH --partition=ckpt
#SBATCH --requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=4500M     # measured MaxRSS 2.26 GB (job 42125708, 50 pairs) x 2
#SBATCH --time=00:20:00  # measured 8:25 for 50 pairs + buffer
#SBATCH --output=logs/uq-s2-%A_%a.out
#SBATCH --error=logs/uq-s2-%A_%a.err
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export ROOT=$HOME/a-s-omics/uq DRAWS=${DRAWS:-/gscratch/csde/azimmer/uq-draws} OUT=${OUT:-/gscratch/csde/azimmer/uq-out} BATCH=${BATCH:-50}
mkdir -p "$DRAWS" "$OUT"
module load apptainer/1.1.5
apptainer exec --bind /gscratch --env R_LIBS_USER=$HOME/Rlib-ctrlB $HOME/containers/geospatial.sif \
  Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); source(file.path(Sys.getenv("ROOT"), "uq-stage2.R"))'
