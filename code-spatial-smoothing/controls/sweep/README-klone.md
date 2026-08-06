# Running the sweep on klone

Design and costing are in `SCOPE.md`. This file is the runbook, in execution order.

klone login is Duo / keyboard-interactive only, so it cannot be driven non-interactively from
the dev box. Either run the commands below yourself, or open a multiplexed session once and
they can be run for you: the `hyak` entry in `~/.ssh/config` sets `ControlMaster auto` with
`ControlPersist 4h`, so a single `ssh hyak true` leaves a socket that later commands reuse
without re-authenticating.

## 0. Before you connect: check what you are sending

`m-est-min-matrix.rds` is 142 MB on disk and is gitignored. It is regenerated from
`sam_sandbox/HD Embeddings/m.est.minassay.2025-07-25.Robj` in about 11 seconds, so if it is
missing, rebuild it locally rather than copying it back down.

Home on klone is 10 GB and nearly empty, which is ample. Do not use `/gscratch/csde`, which was
about 94% full.

## 1. Stage

```bash
ssh hyak 'mkdir -p ~/a-s-omics/sweep/{logs,results}'

rsync -av \
  code-spatial-smoothing/controls/sweep/sweep-run.R \
  code-spatial-smoothing/controls/sweep/sweep-ckpt.sbatch \
  code-spatial-smoothing/controls/sweep/m-est-min-matrix.rds \
  hyak:a-s-omics/sweep/
```

## 2. R environment, one time

Needs only `Matrix`, `RcppHNSW`, `igraph`, `elsa`, `raster`. No INLA and no Seurat, so this is a
much lighter environment than Control B needed.

```bash
ssh hyak
module avail r 2>&1 | head            # find the actual module name on the current build
module load r/4.4.0                   # then put this exact line into sweep-ckpt.sbatch
Rscript -e 'install.packages(c("RcppHNSW","igraph","elsa","raster"), repos="https://cloud.r-project.org")'
Rscript -e 'for (p in c("Matrix","RcppHNSW","igraph","elsa","raster")) cat(p, as.character(packageVersion(p)), "\n")'
```

`sweep-ckpt.sbatch` currently guesses at the module line and falls through harmlessly if the
guess is wrong, which means a bad guess shows up as an R-not-found error in the first task's
`.err` file rather than as a failed submission. Fix the line before submitting.

## 3. Submit

```bash
cd ~/a-s-omics/sweep
hyakalloc                             # confirm the account name before relying on it
sbatch sweep-ckpt.sbatch
```

102 array tasks, capped at 40 concurrent, on the `ckpt` partition with `--requeue`:

- task 0: ordered arm over the seven design depths
- task 1: reverse arm over the same
- task 2: ordered arm over all 30 published depths, the reference curve
- tasks 3 to 101: 99 random draws, seven depths each

Each task is single-threaded and takes about 2.5 minutes, except task 2 at about 13 minutes.
Total is roughly 4.3 core-hours, so on `ckpt` this should be minutes of wall time. Every task
rewrites its own output file from scratch, so preemption needs no partial-state handling.

Watch it with:

```bash
squeue -u $USER
tail -f logs/sweep-*_3.err            # message() goes to stderr, not stdout
```

## 4. Collect

```bash
rsync -av hyak:a-s-omics/sweep/results/ code-spatial-smoothing/controls/sweep/results/
```

Then locally:

```r
library(data.table)
r <- rbindlist(lapply(list.files("results", "^sweep-.*csv$", full.names = TRUE), fread))
ref <- r[arm == "ordered"]; rnd <- r[arm == "random"]; rev_ <- r[arm == "reverse"]
for (d in sort(unique(ref$depth))) {
  o <- ref[depth == d, elsa_d50]; n <- rnd[depth == d, elsa_d50]
  cat(sprintf("depth %5d  ordered %.4f | random %.4f +/- %.4f | reverse %.4f | p(random<=ordered)=%.3f\n",
      d, o, mean(n), sd(n), rev_[depth == d, elsa_d50],
      (1 + sum(n <= o)) / (1 + length(n))))
}
```

Lower ELSA means more spatially coherent domains, so the question is whether the ordered arm
sits below the random distribution and by how much.

## Reading the result

- **Random sits on top of ordered at every depth.** Depth is what matters, feature identity is
  not doing the work, and the non-reducibility claim is clean.
- **Random is clearly worse.** Specific high-variance features carry the result, and the claim
  has to be reworded from "many fields" to "a set of high-variance fields".
- **Random is better.** Variance ordering is actively suboptimal, which would be surprising and
  would want checking before it is believed.

The reverse arm bounds the identity effect from the other side, so ordered and reverse together
say how much room there is between the best and worst orderings at each depth.

## If you would rather not use klone

The same script runs locally with no changes. On 8 of 16 cores it is roughly 35 to 45 minutes,
extrapolated from single-threaded measurement rather than measured directly:

```bash
cd code-spatial-smoothing/controls/sweep
distrobox enter emacs-r -- bash -c 'ASOMICS_SWEEP_ROOT=. Rscript sweep-run.R 0 101 8'
```

Given the job is 4.3 core-hours, staging 142 MB and installing five packages on klone is
arguably more total effort than just running it here. klone wins if you want the machine free,
or if the draw count later goes up by an order of magnitude.
