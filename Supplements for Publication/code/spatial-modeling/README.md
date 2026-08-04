# `spatial-modeling/`: preprocessing, model fitting, convolution

Stages 1, 3 and 4 of the pipeline. Start from `../README.md`, which gives the full run order
across both code directories; this file covers what is specific to these scripts.

A working copy of this code also lives at `code-spatial-smoothing/raw-counts/` in the
repository root, under different file names. It is the version the paper's results were
actually produced from, and it is useful for checking this archived copy against. It differs
in three known ways: it carries commented-out SLURM array-job scaffolding, its feature loop
was left mid-resume at `f.to.mod[882:length(f.to.mod)]`, and its `02` script has no
`if(run.lik.comp)` gate around the model-comparison block. The archived copy here is the
cleaner of the two.

## Contents

| File | Role |
|---|---|
| `0000-process-raw-data.R` | Stage 1. Raw GEM to binned counts plus mesh. Standalone |
| `make.pred.grid.v2.R` | Helper, sourced by the setup scripts. Not used in the paper's configuration |
| `supp-4-lik-comparison/` | Stages 3 and 4. **The paper's main pipeline** |
| `supp-3-raw-v-norm-joint-feat-total-model/` | Side experiment: the joint latent total-count model |

## Run order

```r
setwd("/path/to/a-s-omics")

## Stage 1, once. Set run.entire.preprocess <- TRUE inside on a first run.
source("Supplements for Publication/code/spatial-modeling/0000-process-raw-data.R")

## Stage 3. Edit repo.fp at the top first. Sources 01 for setup, then 02 per feature.
source("Supplements for Publication/code/spatial-modeling/supp-4-lik-comparison/00-launch-experiments.R")

## Stage 4. Needs i.d, o.d and the packages already set; o.d must point at the
## stage-3 output directory. Its header carries the assignments as comments.
source("Supplements for Publication/code/spatial-modeling/supp-4-lik-comparison/03-array-post-process.R")
```

Stage 1 writes into `data-inputs/mouse-embryo-raw/`. Stages 3 and 4 write into
`data-outputs/mouse-embryo-raw/<Sys.Date()>/`, so a rerun lands in a new dated directory
rather than overwriting the archived 2025-04-10 run. Stage 4 must be pointed at whichever
dated directory you want to post-process; it does not infer it.

## Paths you must edit

| File | Line | What to set |
|---|---|---|
| `0000-process-raw-data.R` | `setwd("~/Dropbox/genetics/a-s-omics/data-inputs/mouse-embryo-raw/")` | The input directory on your machine |
| `supp-4-lik-comparison/00-launch-experiments.R` | `repo.fp <- "path/to/repo"` | The repository root |
| `supp-3-.../01-array-run-job.R` | `repo.fp <- "path/to/repo"`, immediately overwritten by `repo.fp <- "~/Dropbox/genetics/a-s-omics/"` | The repository root. Delete or edit the second assignment, which is the one that takes effect |
| `supp-4-lik-comparison/03-array-post-process.R` | commented header block | Uncomment and set `setwd`, `i.d`, `o.d` |

`i.d`, `o.d` and `c.d` in the two setup scripts resolve from `getwd()` and need no editing
once the repository root is set. `c.d` was repaired on 2026-08-04; see the change note at the
bottom of `../README.md`.

## `supp-4-lik-comparison/`: the main pipeline

Four scripts. `00` launches, `01` is setup, `02` is the model, `03` is post-processing.

**Parameters set in `01-load-pkg-set-io-set-params.R`:**

- `agg.size <- 50`, the bin size. The other binned objects from stage 1 exist but are unused.
- `sub.bound`, the square field of view: x from 15,000 to 25,000, y from 11,500 to 18,000.
  At bin 50 that is a 200 x 130 grid, 26,000 bins.
- the feature list, loaded from `features.to.use.2025-04-04.Robj` (918 features, stage 2's
  output, tracked here so stage 3 runs without rerunning stage 2), less
  `f.to.drop <- c("a", "Calcb", "Cckar")`. **915 attempted; 912 produced output.**
- `pred.grid.res <- 200` is set and, as its own comment says, not used.

Three options for choosing features are present in the file. Options 1 and 3 are commented
out; **option 2, the stage-2 feature list, is the live one** and is what the paper used.

**Mesh.** Built in `00-launch-experiments.R` from `sub.bound`, not from the data locations,
so the coarse region begins immediately outside the observed bins. `max.edge` is
`c(0.15, 2) * (x-range / 15)` at bin 50, with `cutoff = max.edge / 10` and
`min.angle = 21`. The inner multiplier is looked up from a small table keyed on `agg.size`.

**Priors, set in `02-array-run-model.R`.** Penalised-complexity Matérn priors,
`inla.spde2.pcmatern`, `alpha = 2`, with a sum-to-zero constraint. Each is given as
`c(a, b, c, d)` meaning `P(range < a) = b` and `P(sigma > c) = d`:

| Component | Prior | Used by |
|---|---|---|
| feature count field | `c(500, 0.95, 0.1, 0.05)` | Poisson, ZIP, ZAP |
| feature presence field | `c(1000, 0.95, 10, 0.05)` | ZAP only |

Note that the count-field prior is set twice: once near the top of the script at
`c(500, 0.95, 1, 0.05)`, then reassigned inside each model block to
`c(500, 0.95, 0.1, 0.05)`. **The second assignment is the one in force**; the first is dead.
Read the values from inside the model blocks, not from the header.

**The three likelihoods.** All three put the observed per-bin total count in as the Poisson
offset `E`, so the linear predictor is on the log rate-per-total-count scale.

- **Poisson**: `feat.count ~ feat.count.int + feat.count.field`, family `poisson`.
- **ZIP**: same formula, family `zeroinflatedpoisson1`, adding INLA's zero-probability
  hyperparameter.
- **ZAP**: two likelihoods fitted jointly. A `binomial` presence model on all bins with its
  own SPDE field, and an `nzpoisson` (zero-truncated Poisson) count model on the non-zero bins
  only.

`bru_max_iter = 1` throughout, so `bru()` performs a single linearisation rather than
iterating to convergence. `n.samp <- 1000` posterior samples for prediction. CPO is enabled
via `bru_options_set` so the PIT can be computed.

Fitting and prediction are each wrapped in `try(..., silent = TRUE)`, and every downstream
block guards on `exists()`. This is how the 3 non-converging features fail without killing the
loop, and it is why their absence shows up as missing output files rather than as an error.

**Scores.** MAE, RMSE, MDS (Dawid-Sebastiani) and MLG (mean log score) per model per feature,
into `{feat}-scores.csv`. Stage 4 averages these across features, excluding rows with
`MAE >= 100` or `abs(MDS) >= 100` as non-convergent; the script's inline comments record that
this excludes 4 and 3 rows respectively. That average is Table 1.

**Convolution, in `03-array-post-process.R`.** For each ligand-receptor pair in
`fantom-hierarchy.Robj` whose ligand and receptor are both among the fitted features, six
operators are computed on four different pairs of inputs, giving the 24 `m.*` columns:

| Operator | Definition |
|---|---|
| `min` | element-wise minimum of the two fields |
| `prod` | element-wise product |
| `gm` | geometric mean, `sqrt(l * r)` |
| `k.01`, `k.5`, `k1` | receptor-ligand equilibrium complex concentration at dissociation constant 0.01, 0.5, 1.0 |

Input variants: `m.data.*` from the raw counts, `m.data.d.*` from the raw densities,
`m.est.*` from the posterior median estimated counts, `m.est.d.*` from the estimated
densities. **The paper uses `m.est.min`.**

Uses the posterior **median** (`l.pred$expect$median`), not the mean.

**Read the pairings from `fantom-hierarchy.Robj`, not from `NICHES::ncomms8866_mouse`.** The
two disagree and only the former gives the paper's 1,477. See `../README.md`.

## `supp-3-raw-v-norm-joint-feat-total-model/`

Not a pipeline stage. This fits the joint model in which the per-bin total count is itself a
latent Poisson field with its own SPDE component, rather than being conditioned on as a fixed
offset. It is the comparison referred to in the manuscript's *The status of the latent
total-count model*.

Same script layout, one step off in numbering: `00` is setup, `01` launches, `02` is the
model. Three differences from `supp-4`:

- Features come from the hand-coded **option 1** list, 13 ligands plus 10 receptors,
  23 features, not from the stage-2 file.
- Both Matérn priors are `c(500, 0.5, 3, 0.01)`, weaker than `supp-4`'s.
- `control.inla = list(int.strategy = "eb")`, empirical Bayes, so the hyperparameter posterior
  is not integrated over.

The manuscript reports this comparison qualitatively ("no visible change, not quantified"),
because the archived S3 File holds joint-model outputs rather than a joint-versus-offset
difference. If a like-for-like comparison is wanted, it has to be run rather than recovered.

## Untracked work in progress

`04-run-diffusion.R` and `diffusion-function-testing.R` in `supp-4-lik-comparison/` explore a
diffusion-based alternative to the convolution operators. They were deliberately not
committed, are not part of the paper, and still carry the aspect-ratio defect corrected
elsewhere. They need `deSolve` and `ReacTran`, which nothing else here does.
