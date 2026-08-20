# Regenerated sub-panels for Figures 3 and 4

For Sam. These are drop-in replacements for the individual maps in the existing Illustrator
files. Each file is the map plus its colour bar, 8 x 5 in at 300 dpi, no baked-in title and
tight margins, so panel labels and layout stay in Illustrator where they already are.

Produced by `export-panels.R`. The merged reference grids `Fig3-raw-vs-normalized.png` and
`Fig4-likelihood-comparison.png` in the parent directory show how they fit together.

## Why they needed regenerating

The published panels were drawn with `asp = qp.res.y / qp.res.x`, which for this field of view
is 130/200 = 0.65. In R base graphics `asp` is the ratio of a y-unit length to an x-unit
length, so 0.65 compresses the y axis against x and **stretches every bin horizontally by
1/0.65 = 1.54**. That is the elongated-pixel artifact Reviewer 2 flagged (R2.6).

The bins are genuinely square: the binning routine steps the x and y breaks by the same
`agg.size`, and the Stereo-seq spot pitch is isotropic. So `asp = 1` is correct and no
caveat to the reviewer is needed. Fixed at source at all 151 occurrences in the pipeline
scripts, so anything regenerated from now on is right.

## Figure 3, raw vs count-normalized topology

Layout is unchanged from the published figure: columns are feature counts, feature counts
divided by total, and total counts; rows are observed, modeled, residual.

| Panel | File |
|---|---|
| r1c1 observed Cd44 counts | `Fig3_r1c1_observed-counts.png` |
| r1c2 observed Cd44 / total | `Fig3_r1c2_observed-normalized.png` |
| r1c3 total counts (nUMI) | `Fig3_r1c3_total-counts-nUMI.png` |
| r2c1 modeled Cd44 counts | `Fig3_r2c1_modeled-counts.png` |
| r2c2 modeled Cd44 / total | `Fig3_r2c2_modeled-normalized.png` |
| r2c3 total counts (nUMI) | `Fig3_r2c3_total-counts-nUMI.png` (same map as r1c3) |
| r3c1 residual, counts | `Fig3_r3c1_residual-counts.png` |
| r3c2 residual, normalized | `Fig3_r3c2_residual-normalized.png` |
| r3c3 | intentionally empty: the total counts are an observed offset, not modeled, so there is no residual |

## Figure 4, likelihood comparison

Rows are Poisson, ZIP, ZAP, raw data; columns are P(non-zero), density (counts/nUMI), feature
counts, residual. The bottom-right panel is the nUMI, as in the published caption.

| Panel | File |
|---|---|
| r1 Poisson, c1 P(non-zero) | `Fig4_r1_Poisson_c1_prob-nonzero.png` |
| r1 Poisson, c2 density | `Fig4_r1_Poisson_c2_density.png` |
| r1 Poisson, c3 counts | `Fig4_r1_Poisson_c3_counts.png` |
| r1 Poisson, c4 residual | `Fig4_r1_Poisson_c4_residual.png` |
| r2 ZIP, c1 P(non-zero) | `Fig4_r2_ZIP_c1_prob-nonzero.png` |
| r2 ZIP, c2 density | `Fig4_r2_ZIP_c2_density.png` |
| r2 ZIP, c3 counts | `Fig4_r2_ZIP_c3_counts.png` |
| r2 ZIP, c4 residual | `Fig4_r2_ZIP_c4_residual.png` |
| r3 ZAP, c1 P(non-zero) | `Fig4_r3_ZAP_c1_prob-nonzero.png` |
| r3 ZAP, c2 density | `Fig4_r3_ZAP_c2_density.png` |
| r3 ZAP, c3 counts | `Fig4_r3_ZAP_c3_counts.png` |
| r3 ZAP, c4 residual | `Fig4_r3_ZAP_c4_residual.png` |
| r4 raw, c1 presence | `Fig4_r4_raw_c1_presence.png` |
| r4 raw, c2 density | `Fig4_r4_raw_c2_density.png` |
| r4 raw, c3 counts | `Fig4_r4_raw_c3_counts.png` |
| r4 raw, c4 total nUMI | `Fig4_r4_raw_c4_total-nUMI.png` |

**Column 1: all 24 panels are now present.** An earlier version of this file said the ZIP and
ZAP "probability of a non-zero value" panels could not be reconstructed without a rerun. That
was wrong, and the correction is worth recording because it rested on checking only the
Poisson object and generalising.

- **Poisson.** P(Y>0) = 1 - exp(-mu), mu the expected count. Exact.
- **ZIP.** There is no spatially varying zero probability to recover. INLA's
  `zeroinflatedpoisson1` carries a *single* zero-probability hyperparameter, so the claim that
  "p_i(s) was never written to disk" mischaracterised the model. The model script stores
  `expect = scaling_prob * lambda * total.count` with `scaling_prob = 1 - p`, so
  `scaling_prob` is recoverable as `expect / (lambda * total.count)`. Checked across all
  26,000 bins it is constant to seven significant figures (coefficient of variation
  9.4e-07), which is what a scalar hyperparameter must look like and confirms the
  derivation. Then P(Y>0) = `scaling_prob * (1 - exp(-lambda * total.count))`.
  For Cd44, `scaling_prob` = 0.999342, i.e. the fitted zero-inflation is about 0.07%.
- **ZAP.** `presence` was in the saved object all along; only the Poisson object is limited to
  `lambda`, `expect` and `obs_prob`. The ZAP presence component is its own SPDE field, it
  varies over space (sd 0.039, range 0.58 to 0.998), and because the count component is
  zero-truncated, P(Y>0) is exactly that presence probability.

No refit was needed, and none was run. Note that the Cd44 P(non-zero) maps are close to
saturated for all three models, because Cd44 is expressed nearly everywhere in this field of
view; that is a property of the feature, not of the reconstruction.

## One more thing worth fixing while you are in there

Figures 5 and 6 have a **separate and unflagged** version of the same problem. Their spatial
panels come from ggplot in `Sup6 - ELSA Analysis.R` and `Sup7 - Clustering vs Features.R`,
and neither script uses `coord_fixed()`. Without it ggplot stretches the panel to fill
whatever device was requested (7 x 5 in for the Sup6 maps, 18 x 6 in for the Sup7 clustering
maps) rather than to the true 200:130 geometry, so those maps are also not at 1:1. Adding
`coord_fixed()` to the spatial maps in both scripts fixes it. Reviewer 2 only flagged
Figures 3 and 4, so this is discretionary, but it is the same defect.

## Still needed from you

- Final composite assembly and panel labelling for Figures 3 and 4.
- The R1.3 pass: larger physical size and reorganisation so compared panels sit adjacent.
- Figure 6's in-figure labels read "1481 Features"; the corrected count is 1477.
- The R1.5 graphical abstract.


## Figure 5 panels (`panels-fig5/`)

Built by `regenerate-fig5.R` (maps), `fig5-input-panels.R` (the shared ligand and receptor
fields), and `fig5-summary-panels.R` (the two summary plots). Assembled by
`compose-fig-5.py`.

Figure 5's aspect defect had a different cause from Figs 3 and 4. Its spatial panels come from
ggplot in `Sup6 - ELSA Analysis.R`, which never calls `coord_fixed()`, so ggplot stretched each
panel to whatever device was requested rather than to the true 200:130 bin geometry. Adding
`coord_fixed(ratio = 1)` fixes it.

**The published layout repeated the ligand and receptor panels in each of the three operator
columns.** Those fields are `l.est` and `r.est` and do not depend on the convolution operator;
all eight were verified byte-identical across the columns. They are now shown once, in a left
gutter beside the mechanism they belong to. The width that recovers goes to the interaction
fields, which grow from 1.02 in to 1.64 in, a 60% increase, which is the substance of R1.3 for
this figure. Their individual colour bars were dropped as well: at gutter width a bar consumed
about 40% of the panel while being unreadable, and these panels are contextual rather than
quantitative reads.

| panel | file |
|---|---|
| shared ligand field, per mechanism | `Fig5_input_{MECH}_ligand.png` |
| shared receptor field, per mechanism | `Fig5_input_{MECH}_receptor.png` |
| interaction field | `Fig5_{prod,gm,min}_{MECH}_interaction.png` |
| resulting domains | `Fig5_{prod,gm,min}_domains.png` |
| entrogram | `Fig5_{prod,gm,min}_entrogram.png` |
| summary, ELSA vs radius | `Fig5_summary_elsa-vs-paramD.png` |
| summary, ELSA vs cluster count | `Fig5_summary_elsa-vs-numclust.png` |

MECH is one of WNT4-FZD6, RSPO1-LGR6, DHH-PTCH1, VEGFA-KDR.

**Height, not width, is what binds this figure.** At full 7.5 in width the composite is taller
than the PLOS ONE limit of 8.75 in, so `compose-fig-5.py` scales it to the box, which yields
5.66 x 8.75 in and leaves 1.84 in of permitted width unused. This is the same constraint that
forced the published version to `0.705\linewidth`. Twenty panels in a portrait stack do not fit
the format, and no amount of sizing fixes that; transposing the figure so operators are rows
and mechanisms are columns would, at the cost of changing what the caption describes.
