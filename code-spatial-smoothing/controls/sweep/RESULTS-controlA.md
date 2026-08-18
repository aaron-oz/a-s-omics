# Control A, second attempt: the pairing does not drive the domains

Run 2026-08-06 on klone `ckpt`, job 38201790, 100 array tasks, 300 job steps all COMPLETED.
99 exclusion-aware scrambles plus the true pairing. Script: `control-A-v2.R`.

Both defects that caused the original withdrawal are fixed. The clustering is single-threaded
and seeded, so single runs are interpretable rather than swamped by run-to-run spread. The
scramble permutes receptors and then repairs, by targeted swaps, every pair that coincides
with a real FANTOM5 mechanism and every self-pair, checked against all 2,557 FANTOM5 pairs
rather than only the 1,477 whose partners are both fitted. **Zero real pairs remained in any
of the 99 nulls.** The receptor multiset is preserved, so the null's marginal field statistics
are identical to the observed.

## Result

| quantity | observed | null (99 scrambles) | range | z | p |
|---|---|---|---|---|---|
| ELSA d=50 | 0.0694 | 0.0724 ± 0.0049 | 0.0577 to 0.0905 | −0.60 | 0.250 |
| ELSA d=100 | 0.0576 | 0.0598 ± 0.0041 | 0.0474 to 0.0748 | −0.53 | 0.240 |
| domains | 21 | 21.05 ± 1.66 | 18 to 26 | −0.03 | 0.640 |
| **distinct fields** | **1449** | **1424.05 ± 5.80** | **1408 to 1436** | **+4.30** | **0.010** |

`p` for the first three is p(null ≤ observed), the one-sided test in the favourable direction,
since lower ELSA means more coherent. For the last row it is p(null ≥ observed).

**On spatial coherence the result is a clean null.** The true pairing sits at z = −0.60 inside
the scramble distribution, and 24 of the 99 scrambles produced more spatially coherent domains
than it did. The domain count is indistinguishable, 21 against a null mean of 21.05. This is
not a case of an effect too small to resolve: the observed value is squarely inside the null,
and the null is now clean and the clustering deterministic, so there is nothing left to blame.

**On field redundancy the result is not null.** The true pairing yields 1,449 numerically
distinct interaction fields where scrambles yield 1,424 on average, and **not one of the 99
scrambles reached the observed value**. Real pairings collapse less under the common minimum
than random ones do. So the pairing does carry structure; it simply is not structure that
shows up in the spatial coherence of the resulting clustering.

## What this means

Taken with the feature-depth sweep, which found that any random few hundred interaction fields
recover the same domains, the picture is consistent and narrow:

> The ~20 tissue domains are a property of the collection of 912 ligand and receptor gene
> expression fields. They do not depend on which fields are used, and they do not depend on
> how the fields are paired.

That is a real finding and a useful one for a literature building pipelines on top of pairing
assumptions, but it is considerably narrower than "domain architecture derived from inferred
spatial ligand-receptor patterns", which is what the manuscript currently implies.

Note that the response letter already commits us to this test. Section R2.8 says the
pair-shuffle control "could now be run to a conclusion under the deterministic configuration
described above, which would remove that ambiguity; we have not done so because it requires
re-running every condition." It has now been run, it did remove the ambiguity, and the answer
is no effect. Having said in the letter that the experiment was available, reporting the
result is not optional.

## The caveat that keeps this from being an overclaim

This is a null **for the common minimum operator applied to raw counts**, which is the
configuration the paper uses. It is not a general statement that ligand-receptor pairing
carries no spatial information, and it should not be written as one.

The reason matters. The common minimum returns whichever partner is limiting, so wherever one
partner is globally dimmer the operator discards the other field entirely along with the
pairing that selected it. That is measurable: for 42.6% of mechanisms the minimum comes from
the same partner at more than 95% of locations, and 28 mechanisms collapse to fields
numerically identical to another mechanism's. The separate effective-rank work in
`../../../NEXT-PAPER-DIRECTIONS.md` finds the pairing signal is invisible on raw counts under
either the minimum or the product (z = −0.51 and −1.49) and decisive after standardising
(z = +9.69 for the minimum). Control A is a fourth, independent measurement pointing the same
way.

So the defensible statement is that this pipeline, on raw counts under the common minimum,
recovers tissue architecture from the gene fields and not from the pairing, and that the
operator choice is the likely reason.

## Replicated in Seurat

Run 2026-08-06, klone job 38598943, 20 tasks all COMPLETED, 19 scrambles at full depth using
the same seeds as the substitute run so the draws correspond. Zero real pairs in any null.

| quantity | Seurat observed | Seurat null (19) | z | substitute z |
|---|---|---|---|---|
| ELSA d=50 | 0.1236 | 0.0925 ± 0.0180 | **+1.73** | −0.60 |
| ELSA d=100 | 0.1147 | 0.0841 ± 0.0173 | +1.77 | −0.53 |
| domains | 19 | 18.47 ± 1.02 | +0.52 | −0.03 |
| **distinct fields** | **1449** | **1424.3 ± 6.3** | **+3.91** | **+4.30** |

**The conclusion replicates: the true pairing does not produce more spatially coherent domains
than scrambled pairings.** Neither implementation finds the effect the control was designed to
detect. The sign of the non-significant deviation differs, which is what noise around a null
looks like: the substitute puts the observed marginally better than the null, Seurat puts it
marginally worse, and in Seurat 17 of 19 scrambles were more coherent than the true pairing.
Note that the Seurat trend, if it firmed up with more scrambles, would say the true pairing is
slightly *worse*, not better.

**The redundancy result replicates cleanly**, z = +3.91 against +4.30. The true pairing yields
1,449 distinct fields where scrambles yield about 1,424, in both implementations.

These two fit together mechanistically. Scrambled pairings collapse more mechanisms onto a
shared limiting partner, which lowers the effective dimensionality of the input, and lower
dimensionality marginally helps cluster coherence. So the pairing does carry structure, it
shows up as reduced redundancy, and if anything it makes the clustering slightly harder rather
than easier.

**A useful validation falls out of the observed arm.** This is a full-depth Seurat evaluation
of the paper's exact configuration, and it gives ELSA(d50) 0.1236 with 19 domains against the
published Sup7 full-depth values of 0.1237 and 20 domains. That agreement to four decimal
places, from an independently written script on a different machine, is a stronger
reproducibility check than the manuscript currently claims anywhere.

## Caveats

- 19 scrambles in Seurat against 99 in the substitute, so the Seurat z values are less precise.
- Control B has still only been run in the substitute, and its null cannot be re-clustered
  without redoing roughly 800 core-hours of model refits, because the permuted fields were
  never saved.
- One clustering seed, deterministic given the input.
- ELSA at d = 50 and d = 100 only.
- The 1,477 mechanisms are clustered as 1,477 columns, including the 28 that collapse to
  duplicates, matching how the published analysis was run.

## Files

`results-ctrlA/ctrlA-observed.csv` and `results-ctrlA/ctrlA-scramble-0NN.csv` for NN = 1..99.
Columns: tag, n_mech, n_distinct, remaining_real_pairs, n_clusters, elsa_d50, elsa_d100,
seconds.
