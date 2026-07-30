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

## 2. The convolution operator carries the biology, and common minimum may be the wrong choice

This is the most surprising thing to come out of the revision, and it cuts against a choice
the current paper makes.

The current paper adopts **common minimum** as the preferred convolution because it gave the
most clusters at the lowest ELSA. But when we asked whether the *true* ligand-receptor pairing
produces more structured fields than a scrambled pairing, using a deterministic
clustering-free statistic (effective rank of the correlation matrix across the 1477 convolved
fields, 8,000 sampled locations):

| operator | true pairing | scrambled | z |
|---|---|---|---|
| common minimum | 38.29 | 38.74 +/- 0.64 | -0.71 |
| **product** | **93.44** | **97.05 +/- 1.89** | **-1.92** (p = 0.077, 25 scrambles) |
| geometric mean | 18.27 | 18.28 +/- 0.09 | -0.09 |

Only the product operator retains any signal of the true pairing. Suggestive, not significant,
and worth running properly.

**Why this is plausible.** The common minimum is a gating operator: it returns the limiting
partner. For **42.6%** of mechanisms the minimum comes from the same partner at more than 95%
of locations (56.8% at more than 90%), and the convolved field correlates with the dimmer
partner's own field at median 0.76. So for a large fraction of mechanisms, common minimum
largely discards one of the two fields, and with it the pairing information. The product is
multiplicative and keeps both everywhere.

**The tension this creates.** Common minimum wins on the criterion the paper used (cluster
count and spatial coherence) but appears to lose on the criterion that actually matters for
the biology (does the true pairing carry more information than a random one). Those are
different questions and the paper conflates them. That is a genuinely interesting result and
possibly the spine of the next paper: **choosing a convolution operator by the spatial
coherence of downstream clusters selects for smoothness, not for signal.**

**What to do.** Redo the operator comparison with the pairing-informativeness criterion as the
primary endpoint, across operators including the kinetic rate-equation prototype, and see
whether the ranking inverts. If it does, the common-minimum recommendation in the current
paper needs revisiting in print.

---

## 3. Why the pair-shuffle control was weak, resolved

Recorded so it is not re-litigated. Two candidate explanations were considered.

1. *Unmatched receptors stand in for true ones because the fields are alike.* **Refuted.**
   Pairwise correlation across the 912 fitted fields is median |r| = 0.061; only 1.3% of pairs
   exceed 0.5. Ligand-ligand is median |r| = 0.068 (2.0% above 0.5), receptor-receptor median
   |r| = 0.058 (0.9% above 0.5). Scrambling substitutes a genuinely different field.
2. *The operator discards the pairing.* **Supported**, see section 2 above.

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
