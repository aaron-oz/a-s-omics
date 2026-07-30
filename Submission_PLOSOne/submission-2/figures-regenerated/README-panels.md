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
| r2 ZIP, c1 P(non-zero) | **not available, see below** |
| r2 ZIP, c2 density | `Fig4_r2_ZIP_c2_density.png` |
| r2 ZIP, c3 counts | `Fig4_r2_ZIP_c3_counts.png` |
| r2 ZIP, c4 residual | `Fig4_r2_ZIP_c4_residual.png` |
| r3 ZAP, c1 P(non-zero) | **not available, see below** |
| r3 ZAP, c2 density | `Fig4_r3_ZAP_c2_density.png` |
| r3 ZAP, c3 counts | `Fig4_r3_ZAP_c3_counts.png` |
| r3 ZAP, c4 residual | `Fig4_r3_ZAP_c4_residual.png` |
| r4 raw, c1 presence | `Fig4_r4_raw_c1_presence.png` |
| r4 raw, c2 density | `Fig4_r4_raw_c2_density.png` |
| r4 raw, c3 counts | `Fig4_r4_raw_c3_counts.png` |
| r4 raw, c4 total nUMI | `Fig4_r4_raw_c4_total-nUMI.png` |

**The two missing panels.** For the Poisson model P(Y>0) = 1 - exp(-mu) with mu the expected
count, so that panel is exact. For ZIP and ZAP the equivalent needs the spatially varying zero
probability p_i(s), and the saved prediction objects contain only `lambda`, `expect` and
`obs_prob`. p_i(s) was never written to disk. Options: rerun those two models saving p_i(s),
or drop column 1 for the ZIP and ZAP rows and say so in the caption. Rather than substitute
something approximate, they have been left out.

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
