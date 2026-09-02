# Independent review of design v2 (2026-08-27)

Verbatim record of the first external review, which prompted design v3. Saved late: an
earlier message said this was already in the repo, and it was not.

The reviewer was given v2, the state-of-play briefing, and the three results documents, and
asked whether the design answers its own question. Verdict was **do not run as written**.

## 1. The stability limit, and the feature-count asymmetry

Stability as a resolution certificate is a known heuristic with known negative results;
it tracks the uniqueness and symmetry of the clustering optimum, not the truth of the
structure, and at large n almost everything looks stable. Usable as a comparative diagnostic,
but "the finest partition it can actually support" is a stronger reading than stability
licenses.

The asymmetry is worse than a worry and invalidates the comparison as written, in the
direction that matters. The perturbation was "drop 20% of each representation's own columns".
The combination arm's columns are redundant by construction: 28 mechanisms collapse to exact
duplicates, and for 42.6% of mechanisms the minimum returns the same partner at over 95% of
locations, so many interaction fields are relabeled single-gene maps. Drop 20% of 1,477 such
columns and the surviving 80% almost surely still contains a copy of every underlying signal;
drop 20% of 912 gene maps and you delete distinct genes outright. The combination arm looks
more stable at high k for mechanical reasons, and the mechanism inflating its stability is the
exact degeneracy the accumulated evidence indicts. Genes appearing in many mechanisms are also
effectively oversampled in the combination arm's resamples.

The bias can run the other way against the rescaled arm: standardizing equalizes variance,
amplifying noise in dim genes, which can degrade stability for reasons about variance
reweighting rather than biology.

Feature resampling may be the wrong perturbation entirely, since it does not touch measurement
noise. Bin resampling has its own trivializer (adjacent bins are near-duplicates under spatial
smoothing). The posterior is the perturbation this pipeline can supply almost for free.

**Structural fix:** perturb upstream of representation construction, on the shared 912 gene
maps, then build every representation from the same perturbed input.

## 2. Agreement measures

- Report homogeneity with completeness, or use adjusted mutual information. Genuine refinement
  has a signature: homogeneity stays high while completeness drops in a structured way.
  Fragmentation loses homogeneity because fragments straddle reference boundaries.
- Add a refinement-tolerant statistic: map each cluster to its majority reference label, merge,
  compute agreement of the merged partition. Cheap, and directly operationalizes "finer but
  faithful".
- The reference has two limits. Its granularity caps everything external. And it is not
  independent histology; it derives from the same Stereo-seq data, which softens its authority
  as an arbiter.
- The held-out genes are the only instrument that can certify a subdivision finer than the
  reference, which is the contested regime. Holding them in reserve is backwards; they should
  be co-primary. But they need a spatially valid null.

## 3. The decision rule

Coherent in form, underspecified in every part that decides the outcome. "Extends the
stability limit" has no estimator, threshold, or uncertainty, and choosing where a smooth curve
falls away after seeing it is a free parameter worth the whole result. Five replicates cannot
support the rule: "differences inside the spread mean no difference detected" converts low
power into a guaranteed verdict, and "no worse" needs an equivalence margin or it passes by
default. "Matched k" is not operationalized. The realistic outcome is a partial order that the
two-arm rule does not compose over. Pre-registration is closer to theater than protection while
the collapse threshold, matching procedure, equivalence margin, replicate count, resolution
grid, and the 1,477-versus-1,449 column choice remain free.

## 4. Most likely confidently wrong conclusion

The design vindicates the common minimum for an artifactual reason: feature resampling within
each arm's own columns makes the min arm look stable to higher k because its column set carries
duplicated copies of the underlying signals; the no-combination arm destabilizes first; the "no
worse at matched k" guard passes vacuously on 5 replicates; and the rule then declares that the
combination extends the stability limit. The paper's central step is rescued by the very
degeneracy that four independent lines of evidence say is its defect.

Runner-up: promoting held-out validation without a spatial null. With 26,000 spatially
autocorrelated bins treated as independent, essentially any contiguous split shows significant
differential expression, validating fragmentation indiscriminately.

## 5. Better design at the same cost

Keep the resolution sweep, curves-not-points, two-stage arms, determinism check first, and the
commitment to report disagreement. Change:

1. Move the perturbation upstream to the shared inputs.
2. Reuse the neighbor graph across resolutions; the expensive step depends only on the feature
   matrix. Verify by timing before banking it.
3. Promote held-out-gene validation to co-primary with a spatial null. The Control A scramble
   partitions and rotated or reflected label maps are free null families.
4. Pre-specify the numbers, not the intentions.
5. Add effective rank as a free diagnostic. If between-arm stability differences track
   effective rank, that is the artifact signature and the stability result is demoted.
