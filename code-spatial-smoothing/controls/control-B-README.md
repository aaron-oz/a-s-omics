# Control B (R2.8): coordinate-shuffle null on klone `ckpt`

Tests whether the emergent tissue domains reflect real spatial structure or are an
artifact of the spatial smoother. Each permutation scrambles the spot-to-location
assignment of the raw counts, re-fits every gene's Poisson field on the spatially
randomized data (paper's exact `bru` model), then runs the same
convolution -> cluster -> ELSA pipeline. Compare the real reference (perm 0) ELSA to
the shuffled-permutation null: if the real domains are much more spatially coherent
(much lower ELSA) than any shuffle, the structure is real, not smoothing.

## Files to stage on klone (into `$HOME/a-s-omics/`, or edit ASOMICS_ROOT)
- `code-spatial-smoothing/controls/control-B-klone.R`   (the task script)
- `data-inputs/mouse-embryo-raw/fantom-hierarchy.Robj`  (L-R pairings, ~26 KB)
- `data-outputs/mouse-embryo-raw/controls/bin50-subregion-raw.rds`  (compact input, 60 MB)
- `control-B-ckpt.sbatch`

Everything else stays node-local or in home; total home footprint < ~200 MB
(input 60 MB + one tiny summary file per permutation).

```bash
# from the dev box, stage (rsync over the Duo-authenticated ssh):
rsync -avR \
  code-spatial-smoothing/controls/control-B-klone.R \
  data-inputs/mouse-embryo-raw/fantom-hierarchy.Robj \
  data-outputs/mouse-embryo-raw/controls/bin50-subregion-raw.rds \
  hyak:a-s-omics/
scp control-B-ckpt.sbatch hyak:a-s-omics/
```

## R environment (one-time, on a klone build/login node)
Need R with: INLA, inlabru, fmesher, splancs, RcppHNSW, igraph, Matrix, elsa,
raster, data.table, parallel. Either `module load` a suitable R and
`install.packages(...)` + INLA from its repo, or use an apptainer image. Put the
right `module load` / `apptainer exec` line into `control-B-ckpt.sbatch`.
(On klone's normal build, Seurat's annoy would also work, but this script avoids
Seurat entirely, so you don't need it.)

## Submit (ckpt = idle-cycle scavenger, the good-citizen choice)
```bash
cd ~/a-s-omics
sbatch control-B-ckpt.sbatch      # array 0-49%20: perm 0 (real) + 49 shuffles, 20 concurrent
```
Confirm the account/partition first with `hyakalloc`; the script defaults to
`--account=ckpt-csde --partition=ckpt`. Each task is ~1 h on 16 cores; with 49
shuffles this is ~13 core-hours/perm x 50 = ~650 core-hours, minutes-to-an-hour of
wall time on ckpt. `--requeue` handles preemption (tasks are idempotent).

## Collect results
```bash
Rscript -e 'library(data.table); \
  r <- rbindlist(lapply(list.files("data-outputs/mouse-embryo-raw/controls", \
       "control-B-perm.*csv", full.names=TRUE), fread)); \
  setorder(r, perm_id); print(r); \
  obs <- r[perm_id==0]; nul <- r[perm_id>0]; \
  cat("real ELSA d50:", obs$elsa_d50, " null mean:", mean(nul$elsa_d50), \
      " p(null<=real):", (1+sum(nul$elsa_d50<=obs$elsa_d50))/(1+nrow(nul)), "\n")'
```

## Tuning
- `NSAMP` (env, default 1000): posterior draws for the median field. 500 is plenty
  stable and ~1.5x faster; the local validation used 500.
- `--array=0-49%20`: raise the shuffle count for a tighter p-value, or lower `%20`
  to be gentler on the cluster.
- If ckpt queue is deep, `--partition=ckpt-all` widens the node pool.
