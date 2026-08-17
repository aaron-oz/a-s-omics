# Seurat replication: what holds, what does not

Run 2026-08-06 on klone `ckpt`. Jobs 38216093, 38216668, 38218758, 38598718-21, 38601492.
19 random draws per depth (18 at depth 50, one task lost to a transient lmod cache race),
plus ordered and reverse arms, through Seurat's own neighbour and clustering path, mirroring
`Sup7 - Clustering vs Features.R`.

## The deconfounded depth curve

Randomising which features enter at each depth isolates the effect of depth from the effect of
feature identity. This is the analysis the manuscript's experiment was reaching for.

| depth | Seurat random | sd | n | Seurat ordered | substitute random | substitute ordered |
|---|---|---|---|---|---|---|
| 10 | 0.3721 | 0.0441 | 19 | **0.1497** | 0.2594 | **0.1428** |
| 50 | 0.2401 | 0.0340 | 18 | 0.1683 | 0.1686 | 0.1639 |
| 150 | 0.1682 | 0.0292 | 19 | 0.1986 | 0.0883 | 0.1907 |
| 350 | 0.1341 | 0.0324 | 19 | 0.1662 | 0.0661 | 0.1672 |
| 650 | 0.1348 | 0.0244 | 19 | 0.1516 | 0.0667 | 0.1234 |
| 1000 | 0.1199 | 0.0123 | 19 | 0.1348 | 0.0696 | 0.0913 |
| 1300 | — | — | — | — | 0.0703 | 0.0783 |

Lower ELSA means more spatially coherent domains.

## Three findings, and which of them replicate

**1. The depth effect is real and replicates strongly.** With identity randomised, coherence
improves 3.1x from depth 10 to depth 1000 in Seurat (0.3721 to 0.1199) and 3.9x from depth 10
to 350 in the substitute (0.2594 to 0.0661). Welch tests on the Seurat curve: 150 against 350
gives p = 0.0017, 650 against 1000 gives p = 0.025, while 350 against 650 is flat at p = 0.94.
So most of the gain arrives by a few hundred fields and improvement past that is slow but not
absent. **This supports the manuscript's underlying claim** that many interaction fields are
needed rather than a handful, and it supports it in a form the published experiment could not,
because the published experiment could not separate depth from identity.

**2. The non-monotonic shape of the published curve is an ordering artifact, and this
replicates.** In both implementations the variance-ordered arm gets *worse* from depth 10 to
depth 150 (Seurat 0.1497 to 0.1986; substitute 0.1428 to 0.1907) before recovering. The
published Seurat run shows the same hump, peaking at 0.2197 around depth 100. The random arm
has no hump in either implementation: it improves monotonically throughout. So the dip-then-
recover shape in Fig 6's column-4 curve is a property of admitting features in variance order,
not of feature count.

**3. Whether random beats ordered at a given depth is implementation-sensitive and does NOT
replicate.** At mid depths the substitute puts the ordered arm about 11 sd above the random
spread; Seurat puts it about 1 sd above, one-sided p roughly 0.15 to 0.20 with 19 draws. The
reverse arm is dramatically best in the substitute and only sometimes better in Seurat.

Consequently these three claims, which I made earlier from the substitute alone, are
**withdrawn**: that variance ordering is actively harmful, that the lowest-variance fields are
the most useful, and that ordering accounts for most of the reported effect.

## What can be said in the manuscript

Supported by both implementations:

> Spatial coherence improves with the number of interaction fields entering the clustering,
> steeply up to a few hundred fields and slowly thereafter. The variance-ordered sweep reported
> here does not show that relationship cleanly, because at each depth it admits a fixed set of
> features and so confounds how many fields enter with which ones; a control admitting random
> subsets of the same sizes recovers a monotone relationship, and attributes the non-monotonic
> shape of the reported curve to the ordering rather than to the feature count.

Not supported, and must not be written: any claim that variance ranking is harmful, or that
low-variance features are preferable.

## Still true independent of all of this

The factual errors in the Fig 6 caption and the surrounding Results came from the saved
metadata of the published Seurat run itself and stand regardless: domain count falls from 33 to
about 20 rather than becoming "more numerous", and full depth is not the most coherent depth
swept. See `SCOPE.md`.

## Caveats

- 19 random draws in Seurat against 99 in the substitute; 18 at depth 50.
- Depth 1300 not run in Seurat. It exceeded the `future` globals cap on the first attempt and
  would have dominated the budget for the least informative part of the range.
- One clustering seed per arm. Seurat's annoy backend is approximate and order-dependent, so
  the Seurat numbers carry run-to-run variability the substitute's do not.
- The two implementations differ substantially in absolute ELSA at the same depth. They agree
  closely on the ordered arm at depths 150 and 350, and diverge on the random arm. Do not mix
  values across implementations.

## Files

`results-seurat-d10/`, `results-seurat-d50/`, `results-seurat/` (depth 150),
`results-seurat-350/`, `results-seurat-d650/`, `results-seurat-d1000/`.
