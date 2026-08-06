# Scoping: a random-order feature-depth sweep

Written 2026-08-06. Purpose: cost and design a control that separates *how many* interaction
fields enter the clustering from *which* ones, because the published sweep confounds the two.

## Why

The published feature-depth sweep (`Sup7 - Clustering vs Features.R`, Fig 6) enters features in
order of decreasing standardised variance. At depth N you therefore always have exactly the same
N features, the top N by variance. Any depth effect is confounded with the identity of the
features that depth admits.

The control is to draw random subsets of each size and compare. Three arms:

| arm | feature order | what it isolates |
|---|---|---|
| ordered | decreasing variance (as published) | the published curve, re-run in this implementation |
| random | uniform random subsets, R draws per depth | the depth effect with identity averaged out |
| reverse | increasing variance | the opposite extreme, bounding the identity effect |

If the random distribution at depth N sits on the ordered curve, depth is what matters and the
non-reducibility claim is clean. If random is markedly worse, specific high-variance features are
carrying the result and the claim needs rewording.

## Verified timings

Measured on this workstation (16 cores, 60 GB), single-threaded, seeded, on the archived
`m.est.min` fields. One evaluation is one clustering plus one ELSA at d=50. From
`timing-variance-ordered.csv`:

| depth | cluster | ELSA | domains | ELSA(d50) |
|---|---|---|---|---|
| 10 | 5.1 s | 0.0 s | 38 | 0.1428 |
| 150 | 10.3 s | 0.1 s | 45 | 0.1907 |
| 650 | 27.3 s | 0.0 s | 19 | 0.1234 |
| 1300 | 44.0 s | 0.0 s | 23 | 0.0783 |
| 1449 | 48.0 s | 0.0 s | 19 | 0.0761 |

ELSA is free. Clustering dominates and scales close to linearly in depth, about
4.8 s + 0.030 s per feature.

**Implementation validated at full depth.** 19 domains at ELSA(d50) 0.0761 here, against Control
B's published real reference of 20 domains at 0.075. That is the check that this clustering
substitute behaves like the one the manuscript's controls used.

## Cost

Depths: 10, 50, 150, 350, 650, 1000, 1300. Full depth is excluded because a random subset of the
whole set is the whole set. Using measured values where available and interpolating otherwise,
one pass over the seven depths is about 145 s.

| arms | passes | single-threaded | 8 workers (extrapolated) |
|---|---|---|---|
| ordered + reverse | 2 | 4.8 min | under 2 min |
| random, R = 49 | 49 | 2.0 h | ~15 min |
| **random, R = 99** | **99** | **4.0 h** | **~35 min** |
| ordered at all 31 published depths | 31 | 13 min | ~3 min |

**Recommended: R = 99, plus all three arms and the 31-depth reference curve. About 4.3 hours
single-threaded, roughly 35 to 45 minutes on 8 workers.** The 8-worker figure is extrapolated
from single-threaded measurement, not measured; treat it as approximate. 8 of 16 cores leaves the
machine usable.

R = 99 rather than 49 because it puts the minimum attainable permutation p at 0.01 rather than
0.02. The convolution-operator work was left with a p floored at 1/9, and that is worth not
repeating when the extra cost is two hours of unattended compute.

Memory: the field matrix is 26,000 x 1,449 doubles, about 300 MB, and `scale()` doubles it
transiently. Budget 1.5 GB per worker, so about 12 GB for 8. There is 41 GB free.

## Validation run, and an early signal

Tasks 0 (ordered) and 3 (random draw 1) were run locally end to end before the harness was
handed over. Both code paths work, and the ordered arm reproduces the standalone timing run
exactly (38 domains at ELSA 0.1428 for depth 10 in both), confirming the clustering is
deterministic under a fixed seed at one thread.

Measured cost of one full pass over the seven design depths: 151 s, against the 145 s
estimated above. So R = 99 is about 4.2 h single-threaded rather than 4.0.

| depth | ordered ELSA(d50) | random draw 1 |
|---|---|---|
| 10 | 0.1428 | 0.2597 |
| 50 | 0.1639 | 0.1923 |
| 150 | 0.1907 | **0.0949** |
| 350 | 0.1672 | **0.0713** |
| 650 | 0.1234 | **0.0761** |
| 1000 | 0.0913 | **0.0698** |
| 1300 | 0.0783 | **0.0713** |

**This is a single random draw and must not be read as a result.** But from depth 150 upward
that one draw is more spatially coherent than the variance-ordered selection at every depth,
by a wide margin in the middle of the range. If that survives 99 draws it lands in the third
branch below, the one flagged as surprising: variance ordering would be actively suboptimal
beyond the first hundred or so features, and the shape of the published sweep curve would be
substantially a property of the ordering rather than of depth.

That would not damage the non-reducibility claim. It arguably strengthens it, since it would
say almost any few hundred interaction fields recover the architecture and there is nothing
special about the high-variance ones. But it would mean Fig 6's curve is telling a different
story than the one the caption tells, so it is worth knowing before the sweep is written up.

## Two things found while scoping

### 1. There are 1,449 numerically distinct interaction fields, not 1,477

Twenty groups of identical rows, 52 rows in total, so 32 redundant. Four are the known Dropbox
conflicted copies. The other 28 are not a file-handling artifact, they are the common minimum
operator doing what it does:

```
GDF9-ACVR2A | GDF9-BMPR1A | GDF9-BMPR1B | GDF9-BMPR2 | GDF9-FXYD6 | GDF9-ORAI2 | GDF9-TGFBR1
NXPH3-NRXN1 | NXPH3-NRXN2 | NXPH3-NRXN3
DKK1-KREMEN1 | DKK1-LRP5 | DKK1-LRP6
BGN-TLR1 | HSP90B1-TLR1 | VCAN-TLR1
CALM1-ADCY8 | CALM2-ADCY8 | GNAS-ADCY8
```

When one partner is the limiting one at every location, `min(L, R)` returns that partner's field
regardless of who the other partner is, so every mechanism sharing that partner collapses onto the
same field. Seven GDF9 mechanisms are one field. This is the same phenomenon already recorded in
`NEXT-PAPER-DIRECTIONS.md` section 2, where for 42.6% of mechanisms the minimum comes from the
same partner at more than 95% of locations, and it is a directly countable instance of it.

Consequence for the manuscript: "1,477 interaction fields" is right as a count of mechanisms and
wrong as a count of distinct fields. The sweep below uses the 1,449 distinct ones.

### 2. The Fig 6 caption's depth claim is contradicted by the run that produced Fig 6

Read from the saved metadata of the 2025-07-12 m.est.min run, which is the run behind Fig 6.

Domain count **falls** with depth and then plateaus:

```
  10 -> 33     150 -> 24     400 -> 20     650 -> 19    1000 -> 19    1300 -> 20
  50 -> 33     200 -> 25     450 -> 18     700 -> 19    1050 -> 18    1450 -> 19
 100 -> 30     250 -> 22     500 -> 16     750 -> 18    1100 -> 19    1481 -> 20
```

Mean ELSA(d50) is non-monotonic. It starts at 0.1432 at depth 10, gets **worse** to a peak of
0.2197 at depth 100, then improves steadily to about 0.11 by depth 1200 to 1450, and is 0.1237 at
full depth.

The caption says "Domains become more numerous and more spatially organized as features are
added." The first half is backwards: they become *fewer*. The Results text also says clustering in
the full feature space "gave the greatest spatial complexity and coherence of the depths we
swept", and the saved ELSA has depths 1250, 1350, 1400 and 1450 all lower than full depth.

Two caveats before anyone acts on this. The published sweep was unseeded and multithreaded, so
each depth is a single draw from a distribution whose spread was measured at about 0.0045 ELSA
units; the small depth-to-depth wiggles are noise. The large movements are not: 33 to 20 domains,
and 0.22 to 0.11 ELSA, are far outside that.

The honest version of the finding is arguably the better story. Adding interaction fields
consolidates a fragmentary partition into fewer, more spatially coherent territories. That is a
more interesting claim than "more domains", and it is what the data shows.

## Note on comparability

ELSA is not comparable between clustering implementations. At full depth, Seurat gives 0.1237 and
the substitute used here gives 0.0761, with 20 and 19 domains respectively. Both are internally
consistent, and Control B's null comparison is sound because it stays inside one implementation.
But the sweep's output curve cannot be overlaid on Fig 6's column-4 curve, so this belongs as a
new supporting figure with its own axis, not as an addition to Fig 6.

Seurat's own `FindNeighbors` aborts in the `emacs-r` container on a libstdc++ assertion, which is
why the substitute exists. If matching Fig 6 exactly matters more than convenience, that is worth
one attempt to resolve first.

## Files

- `m-est-min-matrix.rds` — mechanism-by-bin matrix plus variance ranks, extracted once from
  `sam_sandbox/HD Embeddings/m.est.minassay.2025-07-25.Robj`. Untracked, 300 MB.
- `time-one-evaluation.R` — the timing run above.
- `timing-variance-ordered.csv` — its output.
