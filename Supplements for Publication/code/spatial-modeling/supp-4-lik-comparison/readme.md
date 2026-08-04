# Stages 3 and 4: fit the three count models, then convolve

This is the paper's main modelling pipeline. It fits a spatial model to each ligand and
receptor feature under three count likelihoods, compares them, and combines the fitted
ligand and receptor fields into ligand-receptor interaction fields.

It supports *Spatial Modeling of Tissues for Morphogenic Field Analysis*
(Osgood-Zimmerman & Raredon), PLOS ONE, PONE-D-26-04349.

## Run it

```r
setwd("/path/to/a-s-omics")
## edit repo.fp at the top of 00-launch-experiments.R first
source("Supplements for Publication/code/spatial-modeling/supp-4-lik-comparison/00-launch-experiments.R")
## then, with i.d and o.d pointing at the stage-3 output directory:
source("Supplements for Publication/code/spatial-modeling/supp-4-lik-comparison/03-array-post-process.R")
```

`00` sources `01` for setup and `02` once per feature. `03` is separate and is run afterwards
against a chosen output directory.

Requires stage 1 (`../0000-process-raw-data.R`) to have produced
`data-inputs/mouse-embryo-raw/pre-fit-obj-binned-50.Rdata`. The feature list from stage 2 is
tracked here as `features.to.use.2025-04-04.Robj`, so stage 2 need not be rerun.

Expect roughly 34 GB of output. Fitting all features is the pipeline's expensive step; it is
embarrassingly parallel across features and the loop in `00` can be sliced across processes.

## Files

| File | Role |
|---|---|
| `00-launch-experiments.R` | Loads data, clips to the field of view, builds the mesh, loops over features |
| `01-load-pkg-set-io-set-params.R` | Input/output paths, packages, bin size, field of view, feature list |
| `02-array-run-model.R` | The model: SPDE fit under Poisson, ZIP and ZAP, prediction, scores, plots, save |
| `03-array-post-process.R` | Collates scores into Table 1; convolves ligand and receptor fields |
| `features.to.use.2025-04-04.Robj` | Stage 2's output, 918 features |

## Full documentation

Parameters, priors, likelihood specifications, the convolution operators, feature counts,
which script makes which figure, and the known limits are in `../README.md` and
`../../README.md`. Do not treat this file as sufficient on its own; in particular, read the
notes there on the ligand-receptor pairing source and on the archived
`convolution-data/` file count before comparing your output to the paper's.
