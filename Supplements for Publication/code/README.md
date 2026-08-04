# Analysis pipeline: how to reproduce the results in the paper

Code for *Spatial Modeling of Tissues for Morphogenic Field Analysis*
(Osgood-Zimmerman & Raredon), PLOS ONE, PONE-D-26-04349.

This file is the entry point. It gives the execution order of every script, the input each
one needs and where to obtain it, the intermediate object each one produces and consumes,
which scripts generate which manuscript figures, and the package versions. Two sibling
files carry the per-directory detail:

- `spatial-modeling/README.md` — the model-fitting and convolution stages (stages 1, 3, 4)
- `omics-analysis/README.md` — the feature-selection, clustering, and figure stages
  (stages 2, 5, 6, 7, 8)

The controls reported in the manuscript's *Controls* section are a separate body of code
with its own documentation, at `code-spatial-smoothing/controls/README.md` in the repository
root. They are not part of the pipeline below; they consume its outputs.

Read the **Honest limits** section at the bottom before planning a reproduction. Several
things in here do not run as shipped, and it is better to know which before starting.

---

## Stage overview

The pipeline is eight scripts in two directories, written by two authors, joined at two
points. Read down the `produces` column and up the `needs` column to see the chain.

| # | Script | Needs | Produces |
|---|---|---|---|
| 1 | `spatial-modeling/0000-process-raw-data.R` | raw GEM (11 GB, public) | `pre-fit-obj-binned-{25,50,75,100}.Rdata` |
| 2 | `omics-analysis/Sup2 - Feature Refining.R` | `pre-fit-obj-binned-50.Rdata` | `features.to.use.2025-04-04.Robj` (918 features) |
| 3 | `spatial-modeling/supp-4-lik-comparison/00-launch-experiments.R` | stages 1 and 2 | `prediction-objects/*-pred-{poi,zip,zap}.rds`, `*-scores.csv` |
| 4 | `spatial-modeling/supp-4-lik-comparison/03-array-post-process.R` | stage 3, `fantom-hierarchy.Robj` | `convolution-data/convolution-calcs-*.csv` |
| 5 | `omics-analysis/Sup5 - HD Embeddings.R` | stage 4 | `assay.list.2025-07-26.Robj` (32 clustered assays) |
| 6 | `omics-analysis/Sup6 - ELSA Analysis.R` | stage 5 | ELSA trends per operator; **Fig 5** panels |
| 7 | `omics-analysis/Sup7 - Clustering vs Features.R` | stage 5 | `object.feature.depth.experiment.2025-07-12.Robj`; **Fig 6** rows 1--5 |
| 8 | `omics-analysis/Sup8 - Morphogen Field Visualization.R` | stage 7 | 1,481 interactive HTML fields; **Fig 6** rows 6--7 |

The two joins between the authors' chains are worth naming explicitly, because they are the
places a reproduction can silently diverge:

- **Stage 2 → stage 3** passes `features.to.use.2025-04-04.Robj`. A copy is tracked in
  `spatial-modeling/supp-4-lik-comparison/`, so stage 3 runs without stage 2 having been rerun.
- **Stage 4 → stage 5** passes a *directory* of CSVs, not a single object. Stage 5 reads
  every file in it. See the 1,481-versus-1,477 note under **Honest limits**.

`spatial-modeling/supp-3-raw-v-norm-joint-feat-total-model/` is a side experiment, not a
pipeline stage. It fits the joint latent-total-count model referred to in the manuscript's
*The status of the latent total-count model*, on 23 hand-picked features rather than the full
set. Run it only if you want that comparison; nothing downstream depends on it.

`spatial-modeling/make.pred.grid.v2.R` is a helper sourced by the stage-3 setup script. It is
not run directly, and in the configuration used for the paper its output is not used: the
model predicts at the observed bin locations, not on a separate prediction grid.

---

## Before you start

### Working directory

Every script resolves its paths from `getwd()` at the repository root. Set it once:

```r
setwd("/path/to/a-s-omics")
```

The scripts as archived contain absolute paths from the machines they were written on
(`~/Dropbox/genetics/a-s-omics/`, `/Users/msbr/Dropbox/a-s-omics/`). Those are listed
individually in the two per-directory READMEs and have to be edited before running. Nothing
detects a wrong path for you; a wrong path surfaces as a file-not-found several minutes in.

### Inputs, and where to obtain them

| Input | Size | Where it comes from |
|---|---|---|
| `data-inputs/mouse-embryo-raw/full-raw-t16-data.tsv` | 11 GB | **Public.** Chen et al. 2022 Stereo-seq mouse embryo (MOSTA), E16.5 section E16.5_E1S3, raw GEM at bin 1. <https://db.cngb.org/stomics/mosta/>. Also mirrored as the manuscript's S1 File |
| `data-inputs/mouse-embryo-raw/ligand-receptor-names.RDs` | small | Tracked in the repository. The union of the ligand and receptor symbols used to subset the raw GEM in stage 1 |
| `data-inputs/mouse-embryo-raw/fantom-hierarchy.Robj` | 26 KB | Tracked in the repository. FANTOM5 ligand-receptor pairings as packaged by NICHES. **Stage 4 must use this file, not `NICHES::ncomms8866_mouse` directly**; see **Honest limits** |
| `spatial-modeling/supp-4-lik-comparison/features.to.use.2025-04-04.Robj` | small | Tracked here. Stage 2's output, so stage 3 can be run without rerunning stage 2 |

Stage 1's intermediates (`partial-raw-t16-data.csv`, `total-raw-t16-data.csv`, the four
`pre-fit-obj-binned-*.Rdata`) are derived, large, and not tracked. They are also archived as
the manuscript's S1 File and under `Supplements for Publication/Raw data object/` and
`Supplements for Publication/Intermediate data objects/` (11 GB and 17 GB respectively).

### The field of view

Stages 3 onward do not model the whole section. They model a square sub-region set in the
stage-3 setup script as

```r
sub.bound <- splancs::as.points(matrix(c(15000, 25000, 25000, 15000, 15000,
                                         11500, 11500, 18000, 18000, 11500), ncol = 2))
```

which at bin 50 is a 200 x 130 grid, 26,000 bins. Every count in the paper downstream of
stage 3 is over these 26,000 bins. Stages 1 and 2 operate on the larger tissue boundary
polygon hard-coded in stage 1.

### Compute

Stage 3 is the expensive one. It fits three models (Poisson, ZIP, ZAP) per feature for 915
features, with 1,000 posterior samples each. The archived run took from 2025-04-10 to
2025-04-18 in wall-clock terms, run in slices rather than continuously, so that span is an
upper bound on the compute rather than a measurement of it. For a single-threaded
per-feature figure, the controls work measured the same Poisson fit at roughly 140 s per
feature (see `code-spatial-smoothing/controls/README.md`); the ZIP and ZAP fits in stage 3 are
additional to that, and were not separately timed.

Stage 7 is documented in its own header as taking about a day: roughly 30 iterations at 1 to
40+ minutes each, growing with feature depth.

Stages 3 and 4 together write about 34 GB. Check free disk before starting.

### Packages

The pipeline attaches, across all stages:

```
data.table  dplyr  glue  ggplot2  Matrix  reshape2  stringr  scales  splancs
fmesher  fields  purrr  INLA  inlabru  Seurat  SeuratObject  elsa  tmap  tidyr
cowplot  doParallel  plotly  htmlwidgets  patchwork  knitr  xtable
BPCells  NICHES  ggthemes
```

INLA is not on CRAN; install from <https://www.r-inla.org/download-install>. NICHES is at
<https://github.com/msraredon/NICHES>. BPCells is at
<https://github.com/bnprks/BPCells>. `raster` is used in stages 6, 7 and 8 without an
explicit `require()`; it arrives attached as a dependency of `elsa`. Likewise `viridis`,
`magma`, `cividis` and `turbo` in stages 3 and 4 arrive via `viridisLite`, a dependency of
`fields`. Neither is a bug, but both mean a partial install will fail later than you expect.

A full `sessionInfo()` for a working environment is in
`environment/session-info-2026-08-04.txt`. Versions verified in that environment on
2026-08-04:

```
R 4.5.2        data.table 1.18.0   INLA 25.10.19    inlabru 2.13.0    fmesher 0.6.1
splancs 2.1.45 Seurat 5.5.1        SeuratObject 5.4.0  ggplot2 4.0.1  dplyr 1.1.4
fields 17.1    elsa 1.1.28         raster 3.6.32    Matrix 1.7.4      tmap 4.2
tidyr 1.3.2    cowplot 1.2.0       doParallel 1.0.17 plotly 4.12.0    glue 1.8.0
reshape2 1.4.5 purrr 1.2.1         stringr 1.6.0     scales 1.4.0     patchwork 1.3.2
knitr 1.51     xtable 1.8.4        htmlwidgets 1.6.4
```

**This is a reference environment, not the environment the paper's results were computed
in.** See **Honest limits**.

---

## Stage by stage

### Stage 1. Raw counts to binned counts

`spatial-modeling/0000-process-raw-data.R`

Reads the 11 GB raw GEM, subsets to ligand and receptor genes, clips to a hand-drawn tissue
boundary polygon (32 vertices, found visually with the Shiny tool preserved at the bottom of
the script), then for each of four bin sizes aggregates counts, fills in explicit zeros at
bins where nothing was measured, joins the per-bin total count, and builds the `fmesher`
mesh used by the SPDE model.

Set `run.entire.preprocess <- TRUE` on a first run. It defaults to `FALSE`, which skips the
blocks that write `partial-raw-t16-data.csv` and `total-raw-t16-data.csv`; the binning loop
below then reads files that do not exist.

Produces, per bin size in {25, 50, 75, 100}: `pre-fit-obj-binned-{n}.Rdata` (holding
`all.dat`, `xy.obs`, `coords`, `domain`, `in.dom`, `mesh`, `boundary.coords`) and
`ligand-receptor-binned-{n}.csv`. **The paper uses bin 50 throughout.**

Everything after line 240 is an exploratory scratch section (picking dense ligand-receptor
pairs by hand, and the commented Shiny boundary tool). It is not part of the pipeline, and it
references an undefined `i.d`, so it errors if the file is sourced top to bottom rather than
run interactively down to line 238. Kept for provenance, since the hand-picked list it
produces is what the stage-3 setup script's "OPTION 1" refers to.

`omics-analysis/Sup1 - Raw-to-Bin.R` is a variant of this stage, over *all* genes rather than
ligands and receptors only, producing the whole-section nUMI figures archived as the S1 File.
It is not on the path to any modelling result.

### Stage 2. Feature selection

`omics-analysis/Sup2 - Feature Refining.R`

Loads `pre-fit-obj-binned-50.Rdata`, restricts to the field of view, reshapes to a
gene-by-bin matrix (973 genes x 26,000 bins), and computes two per-gene summaries: the
percentage of bins with a non-zero count, and the maximum count in any bin. Keeps genes at
or above 0.25% positive bins **and** a maximum count of at least 2.

Produces `features.to.use.2025-04-04.Robj`, **918 features**, and the 973 per-gene diagnostic
PNGs archived as the S2 File.

Also produces `mechs.to.use.2025-04-04.Robj`, 1,511 mechanisms, by subsetting
`NICHES::ncomms8866_mouse` to pairs whose ligand and receptor are both among the 918. This
object is *not* what stage 4 uses, and its count is not the paper's count. See
**Honest limits**.

### Stage 3. Fit the three count models per feature

`spatial-modeling/supp-4-lik-comparison/`, scripts `00`, `01`, `02`. Run `00`; it sources
`01` for setup and `02` once per feature.

- `01-load-pkg-set-io-set-params.R` sets `i.d` (inputs), `o.d` (outputs, dated `Sys.Date()`),
  `c.d` (this directory), attaches packages, sets `agg.size <- 50`, sets `sub.bound`, and
  loads the stage-2 feature list. It then drops three features by name,
  `c("a", "Calcb", "Cckar")`, which had failed to fit: **918 - 3 = 915 attempted.**
- `00-launch-experiments.R` loads the binned object, clips to `sub.bound`, builds the mesh at
  the bin-50 resolution, and loops over features.
- `02-array-run-model.R` is the model itself: a PC-prior Matérn SPDE field with `alpha = 2`
  and a sum-to-zero constraint, fitted by `inlabru::bru()` under three likelihoods, with the
  observed per-bin total count entering as the Poisson offset `E`. It then predicts with
  1,000 posterior samples, computes MAE, RMSE, MDS and MLG plus the PIT, and writes.

**Of the 915 attempted, 912 produced output**, which is the paper's fitted-feature count. The
3 that did not are additional to the 3 dropped a priori.

Per feature, into `data-outputs/mouse-embryo-raw/<date>/`: `{feat}-scores.csv`,
`{feat}-diagnostic-plots.png`, two 4x4 comparison PNGs, and three prediction objects in
`prediction-objects/`. The archived run is dated **2025-04-10** and holds 912 score files and
2,736 prediction objects (912 x 3), 17 GB. The comparison PNGs are archived as the S4 File.

Table 1 in the manuscript is the mean of the per-feature scores across this run, excluding
non-convergent rows by the `MAE < 100 & abs(MDS) < 100` filter applied in stage 4.

### Stage 4. Convolve ligand and receptor fields into interaction fields

`spatial-modeling/supp-4-lik-comparison/03-array-post-process.R`

Not self-contained: it expects `i.d`, `o.d` and the packages to be set already, and `o.d`
must point at the stage-3 output directory you want to post-process. Its header carries the
required assignments as comments; uncomment and edit them, or source it after `01`.

Collates the 912 per-feature score files into the Table 1 summary. Then, for each
ligand-receptor pair in `fantom-hierarchy.Robj` whose ligand and receptor are both among the
fitted features, reads the two Poisson prediction objects and combines their posterior median
fields under six operators: element-wise minimum, product, geometric mean, and a
receptor-ligand kinetic equilibrium at each of three dissociation constants (0.01, 0.5, 1.0).
Each operator is applied four ways: to the raw counts, to the raw densities, to the estimated
counts, and to the estimated densities. That is the 24 `m.*` columns per output file.

**The paper uses `m.est.min`**: the element-wise minimum of the two estimated count fields.

Produces `convolution-data/convolution-calcs-{MECHANISM}.csv`, one per pair, and
`convolution-plots/convolution-comparison-poi-{L}-{R}.png`. The archived 2025-04-10 run holds
**1,477** distinct mechanisms; the data directory contains 1,481 files, and the four extra
are Dropbox conflicted copies rather than mechanisms. See **Honest limits**.

### Stage 5. Cluster each representation

`omics-analysis/Sup5 - HD Embeddings.R`

Reads every CSV in the stage-4 `convolution-data/` directory and pivots it: where stage 4
wrote one file per mechanism with 32 columns, stage 5 builds one matrix per column, each
mechanisms-by-bins. Layers whose values are all zero are dropped per matrix. For the four
ligand-side layers and the four receptor-side layers it collapses the mechanism rows back to
unique ligand or receptor names, since a gene appearing in many mechanisms is otherwise
repeated.

Each of the 32 matrices then goes through the same Seurat pipeline: `ScaleData`,
`FindVariableFeatures`, `RunPCA(npcs = 100)`, `RunUMAP(dims = 1:100)`,
`FindNeighbors(dims = 1:100)`, `FindClusters(resolution = 1)`.

Produces one `{layer}assay.2025-07-25.Robj` per layer plus the combined
`assay.list.2025-07-26.Robj`, and a spatial-plus-UMAP PNG per layer. Archived as the S5 File
(38 GB, the bulk of it the 1,481 CSVs copied alongside).

Note that this stage clusters on the top 100 principal components. The manuscript's headline
clustering, in stage 7, does **not**: it clusters in the full feature space with no PCA. The
two are different analyses and should not be conflated.

### Stage 6. Compare convolution operators by spatial coherence

`omics-analysis/Sup6 - ELSA Analysis.R`

For each of the 32 clustered assays, rasterises the cluster labels onto the bin grid and
computes the entropy-based local indicator of spatial association (ELSA, `elsa` package,
`categorical = TRUE`) at 20 neighbourhood radii, `seq(50, 1000, 50)`. Lower ELSA means more
spatially coherent domains.

Produces `assay.raw.elsa.values.2025-07-27.Robj`, `assay.mean.elsa.trends.2025-07-27.Robj`,
`scatter.dat.2025-08-07.Robj`, the two summary PDFs `elsa.v.paramD.pdf` and
`elsa.v.numclust.pdf`, and per-assay entrogram and distribution PNGs. Archived as the S6 File
(1,330 files).

**This is the source of the manuscript's Fig 5 panels**, and of the result that common
minimum gave the most spatially coherent output of the operators tested. The
`ConvolutionPlot()` function at the end of the script draws the ligand, receptor, and
convolved field triplets.

### Stage 7. Feature-depth sweep

`omics-analysis/Sup7 - Clustering vs Features.R`

Takes the `m.est.min` assay from stage 5, ranks features by standardised variance, then
sweeps feature depth: `c(10, seq(50, N, 50), N)`. At each depth it rebuilds the UMAP
embedding and the nearest-neighbour graph **in the full feature space, `dims = NULL`, no
PCA**, clusters at fixed `resolution = 1`, and records the ELSA of the result at 20 radii.

Produces `object.feature.depth.experiment.2025-07-12.Robj` (carrying one
`{n}featureClusters` metadata column per depth), `meta.data`, `embedding.data` and
`graph.data` companions, `clustering.by.nfeatures.raw.elsa.values.Robj`,
`clustering.by.nfeatures.mean.elsa.trends.Robj`, and per-depth PNGs. Archived as the S7 File.

**This is the source of Fig 6's top five rows.** The published domain map is the
`1481featureClusters` column, i.e. the full-depth clustering.

Reproducibility of this clustering is discussed in the manuscript's *Reproducibility of the
clustering*. In short: the approximate nearest-neighbour graph construction is
order-dependent under multiple threads, so the result is exactly reproducible only once the
thread count **and** the random seed are both fixed. Neither condition suffices alone. The
published run fixed neither. The seed changes the answer, so it is a parameter to report.
This script as archived sets no seed and does not restrict threads.

The `makeCluster(8)` / `registerDoParallel(cl)` block near the top is dead code: the loop
below it is a plain `for`, and the `foreach(...) %dopar%` alternative is commented out. The
sweep runs serially. The cluster is also never stopped.

### Stage 8. Interactive 3D interaction fields

`omics-analysis/Sup8 - Morphogen Field Visualization.R`

Loads the stage-7 object, cuts a 2,000 x 2,000 unit sub-region, and writes one standalone
`plotly` surface per interaction field as a self-contained HTML widget. Archived as the S8
File: 1,481 HTML files, 5.5 GB.

**This is the source of Fig 6's bottom two rows**, the eight selected fields rendered as
density surfaces.

`FindAllMarkers(obj, logfc.threshold = 1)` at line 48 is exploratory; its output is assigned
to `genes.to.study` and then not used, because line 77 overwrites the `goi` selection with
all rownames.

---

## Which script made which figure

| Manuscript item | Made by | Notes |
|---|---|---|
| Fig 1, Fig 2 | no script | Schematics, drawn by hand |
| Fig 3, Fig 4 | stage 3, `02-array-run-model.R` | Per-feature `quilt.plot` panels. Composite assembly was by hand in Illustrator; no layout script exists. Regenerated panels at correct aspect are in `Submission_PLOSOne/submission-2/figures-regenerated/panels/` with a placement table |
| Fig 5 | stage 6 | ggplot panels; composite assembled by hand |
| Fig 6 rows 1--5 | stage 7 | |
| Fig 6 rows 6--7 | stage 8 | |
| Table 1 | stage 4, collation block | Reproduces exactly from the archived per-feature score files |
| S10 Fig | `Submission_PLOSOne/submission-2/figures-regenerated/make-controls-figure.R` | The controls figure; has a generating script, unlike Figs 3--6 |
| S9 Table | `Submission_PLOSOne/submission-2/plos-port/S9-notation-to-code.tex` | Maps each symbol in the model specification to the code that implements it |

The composite assembly of Figs 3 through 6 was done by hand and no layout script survives.
This is stated rather than implied: the panels are reproducible, the composites are not.

---

## Honest limits

Things that will cost a reproducer time if they are not said plainly.

**The environment of the original run was not recorded.** The prediction objects carry no
package-version metadata, the fitted model objects were deliberately not saved (about 2 GB
each), and no `sessionInfo()` was captured at the time. `environment/session-info-2026-08-04.txt`
is a working environment as of 2026-08-04, verified to resolve the pipeline's dependencies and
used for the manuscript's controls, but it is **not** the environment that produced the
2025-04-10 fits. INLA in particular changes its internal defaults between versions, so exact
numerical reproduction of Table 1 from a fresh fit should not be assumed. Table 1 does
reproduce exactly from the archived per-feature score files, which is a weaker claim and the
one we can support.

**Absolute paths from the authors' machines are baked into the scripts.** The stage-3 and
supp-3 setup scripts pointed `c.d` at
`code-spatial-smoothing/raw-counts-for-nature-methods-submission/...`, a directory that does
not exist in this repository, which made the archived code fail at its first `source()`. Every
stage-5 through stage-8 script begins with a `setwd("/Users/msbr/...")`. See **Changes made
on 2026-08-04** below for what was corrected and what was left alone.

**The pairing source matters, and the two available sources disagree.** Stage 4 reads
`fantom-hierarchy.Robj` and yields 1,477 mechanisms among the fitted features. Stage 2's
`mechs.to.use` object reads `NICHES::ncomms8866_mouse` and yields 1,511, over the 918
pre-fitting features. Using `NICHES::ncomms8866_mouse` against the 912 fitted features gives
1,491. **1,477 is the paper's number and `fantom-hierarchy.Robj` is the file that produces
it.** Substituting the NICHES object will look like a discrepancy in your results and is not
one.

**The archived `convolution-data/` directory holds 1,481 files but 1,477 mechanisms.** Four
are Dropbox conflicted copies, named
`convolution-calcs-PDGF{C,D}_PDGFR{A,B} (... conflicted copy 2025-05-15).csv`. Stage 5 reads
every file in the directory, so the archived stage-5 objects and everything downstream of them
carry 1,481 rows, and the published Fig 6 in-image labels read "1481 Features". The correct
count is 1,477, and the manuscript's captions say 1,477. **A fresh run of stage 4 will produce
1,477 files, and the resulting stage-5 objects will differ from the archived ones by those
four duplicated rows.** Delete the conflicted copies before rerunning stage 5 if you want to
match the archived objects; leave them if you want the correct count.

**The clustering in stage 7 is exactly reproducible only under conditions the archived script
does not set.** Fixing the seed alone is not enough, and restricting to one thread alone is
not enough; both are required. The published Fig 6 domains are one realization from a run that
fixed neither, so the specific domain count should be read as approximate. This is described
in the manuscript and quantified in `code-spatial-smoothing/controls/README.md`.

**Stages 5 through 8 have no counterpart in the tracked working tree.** For the
`spatial-modeling` scripts, a working copy also lives at `code-spatial-smoothing/raw-counts/`,
so the archived version can be checked against it. The `omics-analysis` scripts exist only in
this archived form, transcribed from the second author's working directory. They contain
interactive fragments (`View()`, un-assigned plot calls, `readline()` prompts, exploratory
blocks whose results are discarded) and are best read as a record of what was run rather than
as a batch program.

**Two scripts in `supp-4-lik-comparison/` are work in progress and are not part of the
paper**: `04-run-diffusion.R` and `diffusion-function-testing.R`. They are untracked, they
explore a diffusion-based alternative to the convolution operators, and they still carry the
aspect-ratio defect corrected elsewhere. Ignore them.

---

## Changes made on 2026-08-04

Two transcription defects were repaired so that the archived stage-3 and supp-3 code can run.
Both are recorded here so the change is visible and revertible; neither touches the model, the
data, or any reported result.

1. **`c.d` now points at this directory tree.** In
   `supp-4-lik-comparison/01-load-pkg-set-io-set-params.R` and
   `supp-3-.../00-load-pkg-set-io-set-params.R`, `c.d` was set to
   `code-spatial-smoothing/raw-counts-for-nature-methods-submission/<subdir>`, which does not
   exist in this repository. It now points at
   `Supplements for Publication/code/spatial-modeling/<subdir>`, which makes both
   `source(file.path(c.d, "02-array-run-model.R"))` and
   `source(file.path(c.d, '../make.pred.grid.v2.R'))` resolve. The corresponding `source()`
   lines in the two launcher scripts were updated to match.

2. **`run.lik.comp` is now set in the stage-3 launcher.** `supp-4-lik-comparison/02-array-run-model.R`
   gates its model-comparison and plotting block on `if(run.lik.comp)`, but nothing in that
   directory ever defined the flag, so the script errored before reaching its save block. The
   flag is now set in `00-launch-experiments.R` exactly as it already was in the supp-3
   launcher: `TRUE` when all three model flags are `TRUE`, otherwise `FALSE`. The working copy
   at `code-spatial-smoothing/raw-counts/02-array-run-model.R` has no such gate and runs the
   comparison unconditionally, which confirms the gate was introduced during transcription and
   the missing definition was an oversight.

The `/Users/msbr/...` paths in `omics-analysis/` were **not** rewritten. They appear in
`setwd()` calls scattered mid-script rather than in one setup block, several of them switching
directory partway through to control where plots land, and mechanically rewriting them risked
changing which outputs overwrite which. They are inventoried in `omics-analysis/README.md`
instead, so a reproducer can edit them deliberately.
