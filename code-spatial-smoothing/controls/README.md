# R2.8 controls and baselines (PLOS One round-1 revision)

Reviewer 2 asked for controls establishing that the emergent tissue domains reflect
biologically meaningful spatial structure "beyond smoothing, increased feature
dimensionality, or clustering." This directory holds the three controls run in response,
plus one exploratory line of work that was shelved.

All three controls operate on the paper's featured pipeline: a single Day-16 field of
view, binned at 50, 912 features with converged Poisson fits, and the FANTOM5 ligand-receptor
pairings convolved by common minimum.

| Control | Question | Answer | Script |
|---|---|---|---|
| A | Does the true L-R pairing produce more spatially coherent domains than scrambled pairings? | Yes, modestly. Observed ELSA(d50) 0.057 vs a 99-shuffle null averaging 0.066; one-sided p = 0.040 | `control-A-lr-shuffle.R` |
| B | Does the spatial smoother manufacture coherent domains from noise? | No. Refitting under permuted coordinates collapses the domains: real ELSA(d50) 0.075 with 20 domains, versus a 49-permutation null averaging 0.470 with 4 to 12 domains; p = 0.020, the floor for 49 permutations | `control-B-klone.R` (+ `control-B-ckpt.sbatch`) |
| C | Do the emergent domains recover known anatomy? | Yes. ARI 0.401, NMI 0.604, size-weighted purity 77.6%, p = 0.002 | `control-C-concordance.R` |

Lower ELSA means more spatially coherent domains, so in Controls A and B the observed
value should sit below the null distribution.

Controls A and B are complementary rather than redundant. A holds the fitted fields fixed
and scrambles only the pairing, which isolates the pairing's contribution but is a weak
null, since scrambled pairings still reuse genuinely spatially structured gene fields.
B destroys the spatial structure itself and refits from scratch, which is the direct test
of the smoothing concern.

## Run order

Controls A and C reuse the fitted fields from the canonical 2025-04-10 run and need no
model refitting. Control B refits every gene under each permutation and is the only
expensive step.

```
# 0. one-time: assemble the per-gene field matrix from the 2025-04-10 prediction objects
#    (done automatically on the first run of control A; cached thereafter)
Rscript control-A-lr-shuffle.R 99 100          # Control A: 99 shuffles, seed base 100

# 1. Control C needs the MOSTA anatomical annotation for the same section
python extract-annotation.py                   # -> controls/e16-fov-annotation.csv
Rscript control-C-concordance.R                # -> controls/control-C-result.rds

# 2. Control B: extract the compact refit input, then run the permutation array on klone
Rscript control-B-extract-input.R              # -> controls/bin50-subregion-raw.rds (58 MB)
#    stage to the cluster and submit; see control-B-README.md
sbatch control-B-ckpt.sbatch                   # array 0-49: perm 0 = real, 1-49 = shuffles
```

Control B is not runnable on a workstation in reasonable time. One permutation refits 912
Poisson fields at roughly 140 s each single-threaded, about 35 core-hours; the 50-task
array is roughly 1750 core-hours total. It was run on the UW Hyak `klone` cluster's
preemptible `ckpt` partition (job 36956180, all 50 tasks completed 2026-07-09).
`control-B-README.md` has the staging, submission, and collection commands.

The 50 cluster summary files are tracked in `controls/control-B-klone/`, alongside the 50
SLURM stdout logs (`ctrlB-36956180_*.out`) kept as run provenance. They record the job ID,
the compute node, and per-task memory: the array ran across 33 distinct nodes at 16 CPUs
per task, peaking near 2.1 GB of R vector memory per task. The corresponding `.err` files,
which hold the per-permutation refit timings and any convergence warnings, were not
retrieved; `message()` output goes to stderr, so the `.out` files do not contain them.

These files are kept separate from `controls/control-B-perm000.csv`, which is an
independent workstation run of the same real reference permutation at 500 posterior draws
rather than the cluster's 1000. The two agree closely (ELSA(d50) 0.065 locally versus 0.075
on the cluster, 20 domains in both), which is a useful reproduction check but means the
files are not interchangeable.

## Inputs and where they come from

| Input | Size | Provenance |
|---|---|---|
| `data-inputs/mouse-embryo-raw/full-raw-t16-data.tsv` | 11 GB | Public. Chen et al. 2022 Stereo-seq mouse embryo (MOSTA), E16.5 section, raw GEM, <https://db.cngb.org/stomics/mosta/> |
| `data-inputs/mouse-embryo-raw/pre-fit-obj-binned-50.Rdata` | large | Derived from the raw GEM by `Supplements for Publication/code/spatial-modeling/0000-process-raw-data.R` |
| `data-inputs/mouse-embryo-raw/fantom-hierarchy.Robj` | 26 KB | FANTOM5 ligand-receptor pairings as packaged by NICHES. Tracked in this repository so the controls are self-contained |
| `data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/*-pred-poi.rds` | 17 GB | The paper's canonical fitted fields, produced by `Supplements for Publication/code/spatial-modeling/supp-4-lik-comparison/` |
| MOSTA bin-50 annotated `.h5ad` for section E16.5_E1S3 | 5.3 GB | Public, same MOSTA portal. Downloaded to a path outside the repository; `extract-annotation.py` reads it |

## What is tracked here, and what is regenerated

Everything under 1 MB is tracked, so the summary results, domain labels, and figures are
available without rerunning anything. Three large derived objects are gitignored, each
regenerated by a tracked script:

| Not tracked | Size | Regenerated by |
|---|---|---|
| `controls/gene-field-matrix.rds` | 171 MB | `control-A-lr-shuffle.R`, section 1, from the 2025-04-10 prediction objects |
| `controls/bin50-subregion-raw.rds` | 58 MB | `control-B-extract-input.R` |
| `controls/fov-markers-raw.tsv` | 158 MB | `extract-fov-markers.sh` (verified 2026-07-27: regenerates the original byte-for-byte, md5 `6c02e7e8a965ab95cce34db6b5760be6`, 2 min 40 s) |

### Coordinate alignment (Control C)

MOSTA bin names are `"i_j"` bin-50 indices. Raw coordinates are `x = i*50`, `y = j*50`,
with `j` spanning the full section height, confirmed against the `spatial` obsm entry,
which stores `(-j, i)`. Our grid bins are `((x-25)/50, (y-25)/50)`. Under this mapping all
26,000 modeled bins matched an annotation.

## Shelved: the representation-ladder / CytoSignal comparison

`klone-seurat-ladder.R` and the results `representation-ladder.csv`,
`gate-seurat-ladder.csv`, `rescore-purity.csv`, `mosta-plus-annotation.csv`, and
`heart-chamber-resolution.png` belong to an exploratory line of work asking whether the
signalling convolution adds resolution over clustering the ligand and receptor fields
directly. It was shelved on 2026-07-27, for two reasons: no reviewer asked for it, and it
never cleared its own validation gate. The gate required reproducing the 20 published
domains exactly with Seurat's clustering; the best run reached ARI 0.43 against the
published labels, so the ladder comparison mixed two different clustering implementations
and is not interpretable as it stands.

The one result from this line that is clean, because it compares two convolution operators
under a single clustering implementation, is that common minimum yields more spatially
coherent domains than the product operator used by CytoSignal: ELSA(d50) 0.064 versus
0.096.

Treat these files as a record of exploratory work, not as results.

## Known gaps

- `control-A-clean.csv` records a stricter variant of Control A in which the shuffled
  pairings are repaired so that no scrambled pair coincides with a real FANTOM5 mechanism.
  It gives a stronger result than the plain shuffle (p = 0.010 versus 0.040). **Its
  generating script did not survive the session that produced it**, and its observed row
  differs from the plain control's observed row (ELSA 0.0496 with 21 clusters, versus
  0.0572 with 20), which indicates the mechanism list also differed in some way we cannot
  now reconstruct. Do not cite these numbers in the manuscript unless the script is
  rewritten deliberately and rerun. The plain `control-A-lr-shuffle.csv` is fully
  reproducible and is what the response letter cites.
- The Control C figures (`control-C-sidebyside.png`, `control-C-heatmap.png`,
  `control-C-sidebyside-paperdomains.png`) and the domain-label objects
  (`observed-domain-labels.rds`, `paper-domain-labels.rds`) were produced ad hoc and have
  no generating script. The objects themselves are tracked, so downstream work is
  reproducible; the plots would need to be regenerated by a new script for a
  publication-quality figure.
- `paper-domain-labels.rds` holds the 20 published Fig-6 domains, extracted from
  `sam_sandbox/HD Embeddings/m.est.minassay.2025-07-25.Robj` (the `1481featureClusters`
  metadata column, not `seurat_clusters`, which is a different unpublished clustering).

## Environment

Controls A, B, and C were run under R 4.5.2 in the `emacs-r` distrobox container with:

```
data.table 1.18.0   INLA 25.10.19   inlabru 2.13.0   fmesher 0.6.1   splancs 2.1.45
RcppHNSW 0.7.0      igraph 2.3.3    Matrix 1.7.4     elsa 1.1.28     raster 3.6.32
```

`extract-annotation.py` runs under a `uv` virtual environment with `anndata`.

Controls A and B cluster with RcppHNSW plus igraph rather than Seurat's neighbour
backend. Seurat 5.5.1's `FindNeighbors` aborts on a libstdc++ assertion in this container
because RcppAnnoy is built with `_GLIBCXX_ASSERTIONS`. The substitute implements the same
recipe (scale features, clip at 10, no PCA, k = 20 KNN in the full feature space, shared-
nearest-neighbour Jaccard pruned at 1/15, Louvain at resolution 1) and is applied
identically to observed and shuffled runs, so the comparisons are internally valid. On a
standard Seurat build the annoy backend works and Seurat can be used directly.
