# Seurat replication of the sweep: the direction holds, the magnitude does not

Run 2026-08-06 on klone `ckpt`, jobs 38216093, 38216668 and 38218758. 19 random draws plus
ordered and reverse arms, at depths 150 and 350, through Seurat's own neighbour and
clustering path, mirroring `Sup7 - Clustering vs Features.R`.

## Side by side

| implementation | depth | ordered | random (n) | reverse | ordered/random | draws beating ordered |
|---|---|---|---|---|---|---|
| Seurat | 150 | 0.1986 | 0.1682 ± 0.0292 (19) | 0.1967 | 1.18x | 16 of 19 |
| Seurat | 350 | 0.1662 | 0.1341 ± 0.0324 (19) | 0.0874 | 1.24x | 15 of 19 |
| substitute | 150 | 0.1907 | 0.0883 ± 0.0114 (99) | 0.0578 | 2.16x | 99 of 99 |
| substitute | 350 | 0.1672 | 0.0661 ± 0.0087 (99) | 0.0420 | 2.53x | 99 of 99 |

**The ordered arms agree almost exactly**, 0.1986 against 0.1907 at depth 150 and 0.1662
against 0.1672 at depth 350. The two pipelines are measuring the same quantity on the same
input, so the divergence below is not a units or a data problem.

**The direction replicates.** Random subsets are more spatially coherent than the
variance-ordered selection in both implementations, at both depths. Seurat: 16 of 19 and 15
of 19.

**The magnitude does not.** Expressed against the spread of the random draws, the ordered arm
sits at about z = +1.0 in Seurat and about z = +11.6 in the substitute. With 19 draws the
Seurat effect is not significant, one-sided p of roughly 0.15 to 0.20, and more draws would
not change that because the effect size itself is small rather than the sample being thin.

**The reverse arm is mixed.** At depth 350 Seurat does reproduce it, 0.0874 against ordered's
0.1662, which is z = −1.44 against the random spread. At depth 150 it does not, 0.1967 against
ordered's 0.1986, essentially identical. In the substitute the reverse arm was dramatically
best at both depths.

## Verdict

The strong claims from the substitute run do **not** survive replication and must not be
written as results:

- "variance ordering is actively harmful"
- "the shape of the published curve is substantially a property of the ordering"
- "the lowest-variance fields are the most useful"

What survives is a caveat rather than a finding: the reported feature-depth curve is measured
under one particular feature ordering, and randomly chosen subsets of the same size tend to be
somewhat more spatially coherent than the variance-ranked ones, by an amount that this
experiment cannot separate from the spread across subsets. That is worth one sentence in the
manuscript acknowledging that depth and identity are confounded in the experiment as run. It
is not worth a figure and it is not a result.

## What this does not touch

The factual errors in the Fig 6 caption and the surrounding Results were established from the
saved metadata of the published Seurat run itself, not from either sweep implementation, and
they stand regardless of anything here. See `SCOPE.md`.

Control A is also unaffected in its logic, though it too has only been run in the substitute.
Its comparison is observed against 99 scrambles within a single implementation on inputs that
differ only in pairing, and its result is a null. Note the direction of the bias this
replication reveals: the substitute *exaggerates* differences between feature configurations
relative to Seurat, by roughly a factor of ten. An implementation that exaggerates differences
returning a null is a conservative result, so if anything Seurat would be expected to return a
null at least as flat. That is an argument, not a measurement, and replicating Control A in
Seurat would settle it: 100 evaluations at full depth, considerably more expensive than these
because cost grows with depth, or about 19 scrambles for a cheaper directional check.

## Caveats

- 19 random draws in Seurat against 99 in the substitute.
- Two depths, 150 and 350. Depth 1300 was dropped after it exceeded the `future` globals cap
  and would have dominated the budget while testing the weakest part of the effect.
- One clustering seed per arm. Seurat's annoy backend is approximate and order-dependent, so
  the Seurat numbers carry run-to-run variability the substitute's do not.

## Files

`results-seurat/` holds depth 150, `results-seurat-350/` holds depth 350.
