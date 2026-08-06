# Random-order feature-depth sweep: results

Run 2026-08-06 on klone `ckpt`, job 38193317, 102 array tasks, 306 job steps all COMPLETED,
zero errors. Design and costing in `SCOPE.md`.

99 random draws per depth, plus the variance-ordered arm as published and a reverse-ordered
arm. Clustering seed fixed at 42 and single-threaded throughout, so between-arm differences
are attributable to the feature subset and not to clustering RNG. 1,449 numerically distinct
interaction fields (see `SCOPE.md` for why not 1,477).

ELSA is negatively oriented: **lower means more spatially coherent domains.**

## The comparison

| depth | ordered | random mean (sd) | reverse | p(random ≤ ordered) |
|---|---|---|---|---|
| 10 | **0.1428** | 0.2594 (0.0376) | 0.3656 | 0.020 |
| 50 | 0.1639 | 0.1686 (0.0224) | **0.0903** | 0.410 |
| 150 | 0.1907 | 0.0883 (0.0114) | **0.0578** | 1.000 |
| 350 | 0.1672 | 0.0661 (0.0087) | **0.0420** | 1.000 |
| 650 | 0.1234 | 0.0667 (0.0060) | **0.0458** | 1.000 |
| 1000 | 0.0913 | 0.0696 (0.0053) | **0.0590** | 1.000 |
| 1300 | 0.0783 | 0.0703 (0.0048) | **0.0613** | 0.940 |

`p = 1.000` means every one of the 99 random draws was more spatially coherent than the
variance-ordered selection at that depth. At depth 1300, 93 of 99 were.

## What it says

**Variance ranking helps only at the very smallest depth, and is actively harmful past it.**
At depth 10 the ordered selection is genuinely good: only 1 of 99 random draws matched it.
By depth 150 every random draw beats it, and the margin is more than a factor of two.

**The lowest-variance fields are the most useful.** The reverse arm, which admits features in
increasing variance, is better than the random mean at every depth from 50 upward and reaches
0.0420 at depth 350. The ordered arm never gets below 0.0761 anywhere. So standardised
variance is not merely uninformative for this purpose, it is anti-correlated with usefulness.

**Coherence saturates near 350 fields and then degrades slightly.** Under random selection the
mean bottoms out at 0.0661 at depth 350 and rises to 0.0703 by 1300. The reverse arm behaves
the same way, minimum 0.0420 at 350. Only the ordered arm is still improving at full depth,
and it is improving because it began in a hole its own ordering dug.

**Domain count converges to about 20 and stays there.** Under random selection the median is
20 domains at every depth from 350 up, range 16 to 24. This is the most robust thing in the
table, and it is arguably good news for the anatomy result: roughly 20 territories are
recovered from almost any few hundred interaction fields.

## Consequences for the manuscript

The full ordered curve in this implementation runs from 0.1428 at depth 10, up to a worst
point of 0.1948 at depth 250, then down to a minimum of 0.0761 at full depth.

So the Results sentence that clustering in the full feature space "gave the greatest spatial
complexity and coherence of the depths we swept" is **true of the experiment as it was run**.
The problem is that it is true only because of the ordering. A random 350 fields is more
coherent (0.0661) than the variance-ordered 1,449 (0.0761), and the best single random draw
anywhere reached 0.0498.

The Fig 6 caption's "Domains become more numerous ... as features are added" remains false in
both implementations; domain count falls, from 33 to about 20 in the published Seurat run and
from 38 to about 20 here.

## The limitation that matters

**This is not Seurat.** It uses the RcppHNSW plus igraph substitute built for the R2.8
controls, because Seurat's `FindNeighbors` aborts on a libstdc++ assertion in the `emacs-r`
container. The substitute was validated at full depth against Control B, 19 domains at ELSA
0.0761 here against the published 20 at 0.075, and it reproduced bit-identically between this
workstation and klone under a different R and igraph. But its absolute ELSA differs
systematically from Seurat's, 0.0761 against 0.1237 at full depth, so the two are not
interchangeable.

Whether the ordering effect replicates under Seurat is **untested**. That is the one thing
that would make this result safe to build on, and it is worth doing before any of it reaches
the manuscript. Seurat 5.5.1 and RcppAnnoy are already installed in `~/Rlib-ctrlB` on klone,
so the environment exists. A reduced design, the ordered arm plus about 19 random draws at
depths 150, 350 and 1300, would settle it; Sup7's own timings were 1 to 40 minutes per
evaluation, so roughly 60 evaluations is hours of compute, minutes of wall time at 40
concurrent on ckpt. I have not measured Seurat timings on klone, so treat that as an estimate.

## Other caveats

- One clustering seed. Deterministic given the subset, so there is no clustering noise, but
  the absolute values are one seed's worth.
- ELSA at d = 50 is quoted throughout; d = 100 is in the CSVs and behaves the same way.
- The reference curve arm covers 29 depths from 10 to 1,449; the three-arm comparison covers
  the seven design depths only.

## Files

`results/sweep-{ordered,reverse,reference}-000.csv`, `results/sweep-random-0NN.csv` for
NN = 1..99. Columns: arm, draw, depth, domains, elsa_d50, elsa_d100, seconds.
