# Side experiment: the joint latent total-count model

**This is not a pipeline stage.** Nothing in the paper's main results depends on it, and
nothing downstream consumes its output. The main pipeline is in `../supp-4-lik-comparison/`.

It supports *Spatial Modeling of Tissues for Morphogenic Field Analysis*
(Osgood-Zimmerman & Raredon), PLOS ONE, PONE-D-26-04349.

## What it does

The main pipeline conditions on the observed per-bin total count, entering it as a fixed
Poisson offset. This experiment instead treats the total count as a latent spatial field in
its own right, fitting it jointly with the feature field: a Poisson likelihood for the total
with its own SPDE component, and a second Poisson likelihood for the feature whose linear
predictor is the sum of the feature and total predictors. It therefore estimates

- the total counts as a smooth field,
- the feature density, meaning counts per total count,
- and the feature counts, as the product of the two.

This is the comparison referred to in the manuscript's *The status of the latent total-count
model*. The manuscript reports it qualitatively, because the archived S3 File holds
joint-model outputs rather than a joint-versus-offset difference. A like-for-like comparison
would have to be run rather than recovered.

## Run it

```r
setwd("/path/to/a-s-omics")
## edit repo.fp at the top of 01-array-run-job.R first: it is assigned twice,
## and the second assignment is the one that takes effect
source("Supplements for Publication/code/spatial-modeling/supp-3-raw-v-norm-joint-feat-total-model/01-array-run-job.R")
```

`01` sources `00` for setup and `02` once per feature. Note the numbering is one step off from
`../supp-4-lik-comparison/`: here `00` is setup and `01` launches.

Requires `data-inputs/mouse-embryo-raw/pre-fit-obj-binned-50.Rdata` from stage 1
(`../0000-process-raw-data.R`).

## How it differs from the main pipeline

| | `supp-4` (main) | `supp-3` (here) |
|---|---|---|
| Total count | fixed offset `E` | latent field, fitted jointly |
| Features | 918 from stage 2, less 3 dropped | 23 hand-picked, 13 ligands and 10 receptors |
| Matérn priors | `c(500, .95, .1, .05)` and `c(1000, .95, 10, .05)` | both `c(500, .5, 3, .01)` |
| Hyperparameters | integrated | empirical Bayes, `int.strategy = "eb"` |
| Likelihoods | Poisson, ZIP, ZAP compared | joint Poisson-Poisson only |

## Full documentation

Run order, inputs, intermediates, package versions, and known limits are in `../README.md`
and `../../README.md`.
