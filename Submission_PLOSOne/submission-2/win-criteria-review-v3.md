# Independent review of design v3 (2026-09-01)

Verdict: **run with modifications, not as written.** Items 1 to 3 blocking.

The reviewer accepted the surviving architecture (resolution sweep, curves not points,
upstream perturbation, effective-rank guard, staging, the measured blocking checks) and
concentrated on what the two prior rounds had not touched.

## 1. The primary measure is not yet an instrument

**1a. The sibling test decides the bias direction, and it was never defined.** Split a real
territory A into A1 and A2. If "supported" means differing from siblings *pooled*, or from
*any* sibling, both fragments pass and the count rises from 1 to 2 under pure fragmentation,
which is the exact inverse of the property claimed for it. If it means distinguishable from
*every* sibling, the property holds but fineness is penalized hard. So the instrument's
direction is an unspecified analyst choice. Recommended: sibling = spatially adjacent domain,
supported = distinguishable from each adjacent neighbor. Note "sibling" cannot mean
merge-tree sibling, because Leiden partitions at different resolutions are not nested.

**1b. The 99-shift null and Benjamini-Hochberg together build a cliff against fine
partitions.** With 99 shifts the smallest attainable p is 0.01. BH at 0.05 across D domains
rejects domains at that floor only if at least D/5 reach it: a 36-domain partition needs 8, an
8-domain partition needs 2. Finer partitions face a discreteness penalty with no biological
content, on top of genuine per-domain power loss. Fix: 999 shifts, and drop BH for a fixed
per-domain alpha, since the count is a comparative ranking instrument rather than an
inferential claim about any single domain.

**1c. Toroidal shifts need a tissue-mask procedure.** The section is an irregular blob in a
bounding box, not a torus. Shifts move labels onto background and leave tissue bins unlabeled,
so null domains are truncated and split, and the null statistic is computed on different
objects than the observed one. Direction of bias unchecked. Specify mask handling and report
sensitivity to one alternative null family.

**1d. Three further undefined choices on the primary outcome:** whether sub-200-bin domains
count as unsupported or are dropped (if dropped, arms are compared over different fractions of
tissue and an arm can park noise for free, so report a coverage companion); which held-out
genes, since most of 27,600 are near-zero at bin 50 and a filter is unavoidable; and what DE
statistic.

**1e. The measure has never been calibrated.** Before it adjudicates anything, run label
surgery on an existing partition: random splits of true domains (the count must not rise),
random merges (must fall), boundary jitter (must fall). No clustering required.

## 2. The two instruments share a confounder

The independence argument fails on redundancy. A lower-rank column set gives smoother
embeddings and coarser structure, which helps stability, and larger smoother domains, which
pass held-out DE more easily. The project's own Control A data shows the mechanism: scrambles
that collapsed more fields got marginally better coherence. Held-out support is downstream of
clustering quality and redundancy affects clustering quality, so both instruments load on it.
The effective-rank demotion must therefore apply to both, not only to stability.

Unspecified and consequential: whether held-out support is computed on the full-input run or
per perturbed replicate. If per replicate, the instruments share the perturbation and
"agreement" manufactures concord rather than protection. Fix to the full-input run.

Decision rule 1 also contradicts section 4: section 4 declines to estimate where the stability
curve falls away, while rule 1 requires "remains stable to a comparable or greater cluster
count", which is that estimate. And no aggregation over the ten-point grid is specified.

## 3. The reversed asymmetry is real, and v3's claim about it is false

"Every arm now loses the same biological information per replicate" is untrue. Combination
arms keep a mechanism only when both partners survive, about 64% of columns, and a surviving
gene whose partners all died contributes nothing to a combination arm while remaining fully
present in the no-combination arm. Losses also arrive in correlated blocks, since dropping one
hub ligand kills all its mechanisms, which depresses between-replicate agreement through a
second channel. Both effects run against the combination arms, the direction the accumulated
evidence already points, which is when a mechanical bias is most dangerous because it reads as
confirmation.

It cannot be equalized within gene deletion: column survival is quadratic in gene survival, so
equal gene loss forces unequal column loss, and equalizing column loss forces the arms to see
different inputs. The perturbation family with no asymmetry is posterior draws, which v3
dropped without recording a reason. Otherwise, split replicates across two survival
intensities as a ranking-invariance check.

## 4. Two gaps in what Stage 1 can conclude

**4a. Missing factorial control.** The arms are none-raw, min-raw, min-rescaled, product-raw.
The rescaled-minimum arm differs from the baseline in two factors at once, and rescaling is
generically helpful. If it wins, the design cannot say whether the combination step or the
rescaling did it, because no-combination-on-rescaled is not an arm. Add it, about +25% cost.
This is the outcome the authors are hoping for, so its uninterpretability matters most.

**4b. The per-field question is never scored.** Everything runs through clustering, yet the
one place the pairing demonstrably carries signal is per-field: distinct-field count at
z = +4.3 and +3.9, and the effective-rank result at z = +9.69. A design that scores only
clustering can conclude "combination adds nothing" while the combination's value sits in
objects it never measured. Add a clustering-free per-field comparison using the Control A
scrambles already on disk.

## 5. Most likely confidently wrong conclusion

Rescaled common minimum is declared the winner, rescuing the paper's central step, through
artifacts. Standardization amplifies noise in dim genes, producing finer partitions; a
pooled-or-any sibling test rewards fragmentation, so that arm posts the highest supported-
domain count; stability does not veto because curves overlap within bands; "no worse at
matched count" passes inside the margin; and the generic benefit of rescaling is credited to
the combination step because the rescaled no-combination control does not exist. Items 1a, 1e
and 4a sever this chain.

## 6. Smaller items

- Matching on achieved cluster count conditions on an outcome; report matched-resolution too.
- 20 replicates give 190 dependent pairs, so a CI treating pairs as independent is
  anticonservative. Jackknife over replicates.
- The matched-count tolerance lacks a tie rule.
- The 2.1 h graph-build figure extrapolates from two points and could be 2x off.
- Single section, single timepoint; scope all conclusions accordingly.
