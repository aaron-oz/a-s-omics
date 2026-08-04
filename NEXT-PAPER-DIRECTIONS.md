# Directions for the next paper

Working notes, started 2026-07-30 during the PLOS One round-1 revision. These are things that
came out of answering the reviewers which are too large for this manuscript but look like the
substance of the next one. Nothing here is started.

---

## 1. A jointly specified multivariate field model

This is AOZ's proposal, recorded properly: give every feature its own spatial field, let the
fields be **differentially correlated through a block covariance structure**, and constrain the
total N(s) to relate to the sum of the feature fields, rather than fitting 912 independent
univariate models and treating the observed total as a fixed offset.

The standard machinery for "block covariance structure that allows different fields to be
differentially correlated" is the **linear model of coregionalization** (LMC): p fields are
expressed as a loading matrix times r << p shared latent spatial factors, so the cross-field
covariance is structured and estimable rather than free. That is the same idea, and it is what
should be written up, not a free p x p cross-covariance, which would be about 10^6 parameters
and is not a serious option.

**Why it is worth doing.** The current model fits each feature independently, so it cannot
borrow strength across features, cannot express that co-regulated features share spatial
structure, and cannot propagate uncertainty from the total into the features. All three are
real limitations and all three are what a reviewer would ask about next.

**What we learned in the revision that bears on it.**

- Fitting the joint feature-plus-total model on 23 features changed the fitted fields very
  little in pattern (correlation 0.986-0.997, median 0.995), somewhat in level (5.1-13.5%
  of the mean fitted count, median 6.6%), and the posterior uncertainty essentially not at
  all (median SD ratio 0.998, median 95% interval width ratio 0.965). So the gain from
  modelling the total jointly, in *this* dataset at *this* depth, is small.
- That is a statement about depth, not about principle. Median 5,740 counts per bin here.
  The joint treatment should matter much more at lower sequencing depth, smaller bins, or
  sparser platforms, which is exactly where a generalisability claim would be tested.

**A structural point that must be handled.** N(s) is *not* the sum of the modelled feature
fields. The 912 curated ligand and receptor genes carry a median of only **4.7%** of total
nUMI (IQR 4.2-5.6%). A model that constrains N(s) to equal the sum of modelled fields must
therefore either model all ~28,500 genes or carry an explicit "everything else" component.

**A cleaner formulation worth considering.** If the total is the sum of the parts, then
conditional on the total the parts are multinomial. Recasting as a **spatial multinomial**
(compositional) model separates composition from total naturally and dissolves the offset
question entirely, rather than answering it. This may be a better framing than a constrained
sum of Poissons.

**Tractability, rough.** r latent factors over the ~10,200-node SPDE mesh is ordinary for
INLA. The loadings are p x r. The hard part is the likelihood: 912 features x 26,000 locations
is about 24M observations. Large but not absurd on a cluster. Estimating r, and identifiability
of the loadings, are the real methodological work.

---

## 2. The convolution operator, and the raw-vs-standardised trade-off

Reworked 2026-08-04 after AOZ caught two errors in the first pass. Numbers below supersede
the earlier ones.

**Two corrections to what was here before.**

1. The scrambled nulls were contaminated. A naive receptor permutation recreates a mean of
   **73.5 genuine FANTOM5 mechanisms out of 1477 (5.0%)**, range 59-94 over 100 permutations.
   Every "scrambled" null was therefore ~5% real pairs. All numbers below use an
   exclusion-aware scramble (targeted swap repair) that leaves exactly zero real pairs.
2. With the null cleaned, **the product result does not hold up**. It was z = -1.92
   (p = 0.077); it is now z = -1.49 (p = 0.222). The earlier claim that "only the product
   retains signal of the true pairing" was over-stated and should not be repeated.

**The statistic.** Effective rank of the correlation matrix across the 1477 convolved fields,
8,000 sampled locations, no clustering. Each convolved field is z-scored across space, so the
correlation matrix compares spatial *pattern* with level and magnitude divided out. Effective
rank = (sum L)^2 / sum L^2 = p^2/||C||_F^2, since a correlation matrix has unit diagonal.
Higher = the fields span more independent spatial modes, i.e. are more distinct from each other.

| convolution input | operator | true pairing | scrambled (n=8) | z |
|---|---|---|---|---|
| raw counts | common minimum | 38.60 | 39.20 +/- 1.17 | -0.51 |
| raw counts | product | 95.48 | 99.26 +/- 2.53 | -1.49 |
| **standardised** | **common minimum** | **12.38** | **11.47 +/- 0.09** | **+9.69** |
| **standardised** | **product** | **288.07** | **273.09 +/- 6.82** | **+2.20** |

**The finding.** On raw counts the statistic cannot tell a true pairing from a random one,
under either operator. On standardised fields it separates them decisively, most strongly
under the common minimum. Note the sign: the true pairing gives a *higher* effective rank,
meaning the true-paired fields are more distinct from one another, where randomly paired
fields collapse toward a more redundant common mode. The permutation p-value is floored at
1/9 with 8 scrambles, so the z-score is the informative summary; the null sd of 0.09 for
standardised/minimum is remarkably tight, which is itself why z is so large.

**Why standardising matters, mechanistically.** On raw counts the common minimum is a gating
operator that returns the limiting partner, and gene magnitudes span orders of magnitude
(Cd44 mean 12.9, Tnfrsf19 mean 0.17). For **42.6%** of mechanisms the minimum comes from the
same partner at more than 95% of locations (56.8% at more than 90%), and the convolved field
correlates with the dimmer partner's own field at median 0.76. So on raw counts the operator
largely discards one of the two fields, and the pairing information with it. Put both fields
on a common scale and the minimum becomes a genuine "both must be high" gate, at which point
the identity of the partner matters.

**The trade-off this exposes, which the current paper does not acknowledge.** The paper argues
for raw counts because they carry absolute binding density, and selects common minimum because
it maximised cluster count at minimum ELSA. But the pairing information appears to be
recoverable only after standardising. So raw-plus-minimum may be buying the absolute-density
interpretation at the cost of the pairing signal. That is a real methodological tension and a
good spine for the next paper: **an operator chosen by the spatial coherence of downstream
clusters is being selected for smoothness, not for whether it preserves the biology.**

**What to do next.** More permutations (8 is enough for z but not for a usable p). Confirm the
direction of the effect is interpretable rather than an artifact of how effective rank responds
to sparsity. Re-run the actual pair-shuffle control (ELSA on clusters) with standardised input,
which is now cheap since the clustering is deterministic. And test the kinetic rate-equation
operator on the same footing.

## 3. Why the pair-shuffle control was weak, resolved

Recorded so it is not re-litigated. Two candidate explanations were considered.

1. *Unmatched receptors stand in for true ones because the fields are alike.* **Refuted.**
   Pairwise correlation across the 912 fitted fields is median |r| = 0.061; only 1.3% of pairs
   exceed 0.5. Ligand-ligand is median |r| = 0.068 (2.0% above 0.5), receptor-receptor median
   |r| = 0.058 (0.9% above 0.5). Scrambling substitutes a genuinely different field.
2. *The operator discards the pairing on raw counts.* **Supported**, see section 2: the
   effect is invisible on raw input and decisive on standardised input.

A third point, from the same work: the published Control A used a naive scramble, so its
null was also ~5% real pairs. It is already withdrawn, but that is a second independent
reason the withdrawal was right.

Also established: the weak result was **not** an artifact of the clustering step. A
deterministic, clustering-free statistic reproduces it. And the clustering is exactly
reproducible once thread count and seed are fixed, so "clustering is an art" is not the
explanation either.

---

## 4. Other platforms

Recorded from the generalisability discussion. We hold much of the mouse brain data and code.
Before running anything, settle: which platforms, chosen to span genuinely different collection
geometries (grid dimensions, gridded versus ungridded, imaging versus spot-based); and what
each additional dataset is meant to demonstrate, since a second dataset showing the same thing
again adds little. Candidates for what to show: that the transfer conditions stated in the
current Discussion actually hold; that the operator result in section 2 is not
platform-specific; that the joint model in section 1 matters more at lower depth.

---

## 5. The visualiser

Deferred until the above is settled, since it should display whatever we end up wanting to
show. Requirements captured so far: click two points to draw a transect and show how a selected
set of features varies along that line **with uncertainty**; select a box or circle and show the
distribution of feature values within it. S8 is already 1,481 standalone interactive HTML
files, so much of the per-mechanism rendering exists; what is missing is the cross-feature and
transect views.
