# `omics-analysis/`: feature selection, clustering, and the figure scripts

Stages 2, 5, 6, 7 and 8 of the pipeline, plus the whole-section summary figures. Start from
`../README.md`, which gives the full run order across both code directories; this file covers
what is specific to these scripts.

**Read these as a record of what was run, not as a batch program.** They are the second
author's working scripts, transcribed as they were. They contain interactive fragments,
un-assigned plot calls whose output goes to the device, exploratory blocks whose results are
computed and then discarded, and absolute paths from the machine they ran on. Nothing here
was written to be sourced top to bottom on another machine, and none of it will be without
edits. The **Defects** section at the bottom lists the places where a literal rerun fails, so
that they are known rather than discovered.

Unlike `../spatial-modeling/`, there is no second working copy of these scripts elsewhere in
the repository to check against. This archived copy is the only one.

## Contents and order

| Script | Stage | Reads | Writes |
|---|---|---|---|
| `Sup1 - Raw-to-Bin.R` | off-path | raw GEM | whole-section nUMI figures, `meta.compiled.MSBR.2025-07-27.Robj` |
| `Sup2 - Feature Refining.R` | 2 | `pre-fit-obj-binned-50.Rdata` | **`features.to.use.2025-04-04.Robj`**, 973 per-gene PNGs |
| `Sup5 - HD Embeddings.R` | 5 | stage-4 `convolution-data/` | **`assay.list.2025-07-26.Robj`**, 32 per-layer Robj and PNG |
| `Sup6 - ELSA Analysis.R` | 6 | `assay.list.2025-07-26.Robj` | ELSA trends per operator, **Fig 5 panels** |
| `Sup7 - Clustering vs Features.R` | 7 | `assay.list.2025-07-26.Robj` | **`object.feature.depth.experiment.2025-07-12.Robj`**, **Fig 6 rows 1--5** |
| `Sup8 - Morphogen Field Visualization.R` | 8 | `object.feature.depth.experiment.2025-07-12.Robj` | 1,481 HTML surfaces, **Fig 6 rows 6--7** |

Sup1 and Sup2 run before the modelling in `../spatial-modeling/`; Sup5 through Sup8 run
after it. Sup5 must precede Sup6 and Sup7, and Sup7 must precede Sup8. Sup6 and Sup7 are
independent of each other.

Sup1 is off the path to any modelling result. It repeats stage 1's binning over *all* genes
rather than ligands and receptors only, purely to produce the whole-section total-count,
ligand-count and receptor-count figures archived as the S1 File.

## Paths you must edit

Every path below is absolute and points at the second author's machine. They are not
concentrated in a setup block: several `setwd()` calls sit mid-script, switching directory
partway through so that different groups of plots land in different places. **Edit them
deliberately, in order, rather than with a global find-and-replace**, or outputs that were
meant for separate directories will overwrite one another.

| Script | Line | Path |
|---|---|---|
| `Sup1` | 20 | `setwd("/Users/msbr/Dropbox/a-s-omics/data-inputs/mouse-embryo-raw/")` |
| `Sup1` | 82 | `setwd(".../mouse-embryo-raw/msbr-2025-04-24")` |
| `Sup1` | 242 | `setwd("~/Dropbox/a-s-omics/sam_sandbox/plotting for figures")` |
| `Sup2` | 12 | `setwd("/Users/msbr/Dropbox/a-s-omics/sam_sandbox")` |
| `Sup2` | 15 | `load(".../data-inputs/mouse-embryo-raw/pre-fit-obj-binned-50.Rdata")` |
| `Sup2` | 26 | `readRDS(".../ligand-receptor-names.Rds")` — see **Defects** |
| `Sup2` | 29--30 | `load("~/Large Files/chen_stereo_mouse_embryo/...")` — outside the repository |
| `Sup5` | 7, 140, 146 | `setwd(".../sam_sandbox/HD Embeddings")` |
| `Sup5` | 30 | `file.path.input <- ".../data-outputs/mouse-embryo-raw/2025-04-10/convolution-data"` |
| `Sup6` | 12 | `setwd(".../sam_sandbox/ELSA/ELSA 2nd run (July 2025)")` |
| `Sup6` | 44 | `load(".../HD Embeddings/assay.list.2025-07-26.Robj")` |
| `Sup7` | 19 | `load(".../HD Embeddings/assay.list.2025-07-26.Robj")` |
| `Sup7` | 129, 198, 247, 378, 497 | `setwd(".../sam_sandbox")` |
| `Sup8` | 10 | `setwd(".../sam_sandbox/morphogen fields")` |
| `Sup8` | 34 | `load(".../clustering vs features/m.est.min/clustering vs features experiment/object.feature.depth.experiment.2025-07-12.Robj")` |

The Sup8 input path shows that Sup7's output was moved by hand after Sup7 wrote it: Sup7
saves into `sam_sandbox/`, Sup8 loads from
`sam_sandbox/clustering vs features/m.est.min/clustering vs features experiment/`. That move
is not scripted.

## Parameters that matter

**Sup2, feature selection.** Two thresholds, both in the script:
`frac.pos.thresh <- 0.25` (percent of bins with a non-zero count) and
`row.max.thresh <- 2` (maximum count in any single bin). Applied to the 973 genes present in
the field of view. Result: **918 features**.

**Sup5, per-representation clustering.** `RunPCA(npcs = 100)`, then UMAP, neighbours and
`FindClusters(resolution = 1)` on `dims = 1:100`. This stage *does* use PCA. It runs over the
32 value columns stage 4 wrote, transposing them from one-file-per-mechanism to
one-matrix-per-representation.

The four ligand-side layers (`l.d`, `l.est`, `l.est.d`, `l`) and four receptor-side layers
(`r.d`, `r.est`, `r.est.d`, `r`) have their rows collapsed to unique gene symbols first,
because a gene participating in many mechanisms would otherwise appear many times. The 24
convolved layers keep one row per mechanism. Layers that sum to zero for a given mechanism
are skipped for that mechanism, so the matrices are not all the same height.

**Sup7, the feature-depth sweep.** This is the manuscript's headline clustering and it is
configured differently from Sup5:

- features ranked by `var.features.rank`, that is standardised variance from Seurat's
  variance-stabilising transform;
- depths `c(10, seq(50, N, 50), N)`, roughly 30 of them;
- `RunUMAP(..., dims = NULL)` and `FindNeighbors(..., dims = NULL)`, meaning **the full
  feature space with no PCA**, which is the point of the experiment;
- `resolution = 1` held fixed at every depth, so that feature count is the only thing
  deliberately varied.

The published Fig 6 domain map is the `1481featureClusters` metadata column on the saved
object, i.e. the full-depth clustering. (On the counting of 1,481 against the correct 1,477,
see `../README.md`.)

**Sup6 and Sup7, ELSA.** `elsa::elsa(r, d = <radius>, categorical = TRUE)` on the rasterised
cluster labels, at 20 radii `seq(50, 1000, 50)` in bin-grid units. ELSA is negatively
oriented: lower means more spatially coherent. Sup6 sweeps over the 32 representations at
fixed feature depth; Sup7 sweeps over feature depth at the fixed `m.est.min` representation.

**Sup8.** Sub-region x from 18,000 to 20,000 and y from 16,000 to 18,000, rendered as
`plotly` surfaces with `aspectratio = list(x = 1, y = 1, z = 0.075)`, saved as self-contained
HTML widgets.

## Clustering reproducibility

Seurat's `FindNeighbors` builds its nearest-neighbour graph by an approximate method that can
run across several threads, and under concurrent insertion the graph depends on the order in
which threads finish. The result is exactly reproducible only when **both** the thread count
and the random seed are fixed; neither condition suffices alone. These scripts set neither.

The manuscript's *Reproducibility of the clustering* reports the measurements: with four
threads and a fixed seed, three replicate runs gave 19, 19 and 20 domains; unseeded
single-threaded runs gave 19, 22 and 22; single-threaded with a fixed seed, five replicates
were bit-identical at 21 domains, adjusted Rand index exactly 1 between every pair. The seed
changes the answer, so it is a parameter of the analysis and should be reported.

The published Fig 6 domains come from the original multithreaded, unseeded configuration and
are therefore one realization. Anyone rerunning Sup7 should expect a somewhat different
partition, and should fix both the thread count and the seed if they want a repeatable one.
`code-spatial-smoothing/controls/README.md` has the full measurement record.

## Defects

Places where the archived scripts do not run as written. None of these affect the published
results, which were produced interactively; they matter to anyone rerunning.

- **`Sup1` line 209 is a garbled line.** It reads
  `fantom5 <- NICHES::ncomms8866_mousx.yfantom5 <- NICHES::ncomms8866_mouse`, a chained
  assignment whose middle target is a mangled name in a package namespace. It does not
  evaluate. The intent is plainly `fantom5 <- NICHES::ncomms8866_mouse`.
- **`Sup2` line 26 has the wrong case for its input file.** It asks for
  `ligand-receptor-names.Rds`; the file is `ligand-receptor-names.RDs`. This works on a
  case-insensitive filesystem such as macOS and fails on Linux. The call is a bare `readRDS`
  whose value is not assigned, so it is a display-only line and can be deleted; the same file
  is read with the correct case in `../spatial-modeling/0000-process-raw-data.R`.
- **`Sup2` lines 29--30 read from `~/Large Files/`**, outside the repository. The block is a
  sanity check that compares gene counts between two sources and assigns nothing used later.
  It can be skipped.
- **`Sup1`'s `run.entire.preprocess` gate is commented out** at lines 28, 45, 84 and 94, so
  the expensive blocks always run. This is deliberate in the archived state, but it means
  the script cannot be re-entered cheaply the way stage 1 can.
- **`Sup7`'s parallel setup is dead code.** `library(doParallel)`, `makeCluster(8)` and
  `registerDoParallel(cl)` are executed, but the loop below is a plain `for`; the
  `foreach(...) %dopar%` alternative is commented out. The sweep runs serially, and the
  cluster is never stopped. The eight workers sit idle for the duration, roughly a day.
- **`Sup7` line 353 plots `to.plot`, not `to.plot.trends`.** At that point `to.plot` still
  holds the per-bin cluster and coordinate frame from the loop above, which has no
  `feat.num`, `elsa.value` or `elsa.dist` column, so `p1` fails when drawn. It is used only in
  `m.est.min feature depth experiment summary trends.png`; the per-depth composites that
  reach Fig 6 rebuild the same panel correctly as `p.a`/`p.a.i` further down.
- **`Sup8` line 48 computes `FindAllMarkers(obj, logfc.threshold = 1)` into `genes.to.study`
  and never uses it.** Line 76 sets a hand-picked `goi` list and line 77 immediately
  overwrites it with all rownames. Both lines are exploratory leftovers; the loop renders
  every field.
- **`raster` is used without being attached explicitly** in Sup6, Sup7 and Sup8
  (`raster()`, `rasterFromXYZ()`). It arrives as a dependency of `elsa`. Not a fault, but an
  install that lacks `elsa` fails in a confusing place.
