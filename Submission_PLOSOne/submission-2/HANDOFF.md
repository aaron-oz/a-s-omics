# Handoff: PLOS One revision, PONE-D-26-04349

Written 2026-08-04. Deadline **2026-08-07**. Everything is committed and pushed at
`0e602fd`; `main` and `plosone-revision` and both remotes are in sync, working tree clean.

**Read this file for context, then work from
`response-to-reviewers.tex` -> the "Open items" section on page 2.** That section is the
authoritative task list and it is checked mechanically. This file exists for what the task
list cannot carry: why things are the way they are, and which of them are fragile.

---

## Where the revision stands

All 23 reviewer items are answered in the manuscript. **Nothing is blocked on the authors.**
Five items carry open markers, all of them production work rather than science:

| item | state | what remains |
|---|---|---|
| JR-5 | INCOMPLETE | S9 Table and S10 Fig are built and staged. S1-S8 are 53.5 GB across 12,275 files and cannot be uploaded; needs an author decision between three options set out in `supporting-information/README-supporting-information.md` |
| R1.3 | TODO | Figures larger and reorganised. Raredon holds the Illustrator files |
| R1.5 | TODO | Graphical abstract. Raredon |
| R2.6 | INCOMPLETE | Cause found and fixed at source; final composite assembly still needs a person |
| R2.9 | INCOMPLETE | Repository documentation for the main pipeline. The reviewer will check this |

Verify with `python3 check-response-letter.py` (exits non-zero on any problem; also prints
which items still hold markers).

## Build

`README-build.md` is accurate and was re-verified. Four artifacts, all currently 0 errors and
0 undefined references: manuscript 28 pp, track-changes 30 pp, S9 Table 4 pp, letter 21 pp.
Three things that will bite otherwise:

- **`latexdiff --math-markup=off` is required.** The documented `--math-markup=whole` produces
  165 errors since the model equation was rewritten for R2.1.
- **Do not count markers with `grep`.** The obvious regex silently matches nothing, because the
  leading backslash-a is consumed as an escape, so a broken gate looks like a passing one. Use
  `check-response-letter.py`.
- **`review-copies/` holds ghostscript-downsampled snapshots** so co-authors can read the
  revision from GitHub; the real build outputs under `plos-port/` stay gitignored. The
  snapshots drift if sources change without refreshing them. They are review copies, **not**
  submission copies.

## Things that are fragile or easy to get wrong

- **Figures 5 and 6 are temporarily scaled down** (`[p]` placement, `\captionsetup{font=small}`,
  width ~0.72\linewidth) so their long new captions fit. This works *against* R1.3 and is a
  stopgap for the reviewable PDF only. Undo it when the figures are regenerated or stripped.
- **Fig 6's in-image labels still read "1481 Features"**; the corrected count is 1477 and the
  caption already says so. Fix when regenerating.
- **The 22 sub-panels for Figs 3 and 4** are in `figures-regenerated/panels/` at correct 1:1
  aspect for Illustrator placement, with a placement table in `README-panels.md`. Two panels are
  deliberately absent: ZIP and ZAP "probability of a non-zero value" cannot be reconstructed,
  because the saved prediction objects hold only `lambda`, `expect` and `obs_prob`, and the
  spatially varying zero probability was never written to disk.
- **Figs 5 and 6 have the same aspect defect, unflagged by the reviewer.** Their panels come
  from ggplot in `Sup6 - ELSA Analysis.R` and `Sup7 - Clustering vs Features.R`, neither of
  which calls `coord_fixed()`. Discretionary, but it is the same bug.
- **The FigShare links are private share URLs, not DOIs.** They can be revoked and would leave
  dead citations. Worth converting regardless of which SI option is chosen.

## Claims that were corrected during the revision, and must not drift back

Each of these was wrong in an earlier draft and is now right. If you find yourself restating
the old version, stop.

1. **The clustering is NOT non-deterministic.** It is exactly reproducible with a single thread
   *and* a fixed seed (five replicate runs, bit-identical, ARI 1.0). Neither condition suffices
   alone. The seed changes the answer, so it is a reportable parameter. Earlier drafts called
   the pipeline non-deterministic, which was a falsehood about the method rather than about the
   configuration used.
2. **Poisson wins on all four scores in Table 1, not three of four.** The MDS percentages that
   suggested otherwise were a ratio of two negative numbers. Those percentages have been removed
   from the table entirely.
3. **`E(F_i(s)) = N(s)lambda_i(s)` is wrong unconditionally**, because lambda is itself a random
   field. It is the conditional mean of the data stage and is written that way.
4. **Treating nUMI as an offset is a conditioning choice, not a claim that N(s) is precise.**
   An earlier draft justified it by the 1.32% Poisson relative sd, which addresses only the
   conditional sampling noise and not the biological variability of the total. The justification
   is now empirical: pattern correlation 0.986-0.997, level difference median 6.6%, and
   posterior uncertainty essentially unchanged (median SD ratio 0.998).
5. **The MOSTA annotations are independent in the sense that matters** (different group,
   different method, our pipeline never saw them). An earlier draft over-caveated this into
   near-uselessness by demanding a second measurement modality, which no concordance check could
   ever satisfy. The paper does state how they were derived, which is spatially constrained
   clustering plus manual labelling, for transparency rather than as an indictment.
6. **Cd44 was chosen for clarity, not because it is representative.** It shows the
   normalisation effect against the developing spine, an unambiguous landmark. The text says so
   and points to S2 and S3 File so a reader can check how general the pattern is.

## Verified numbers, so they are not re-derived

918 curated features -> 3 dropped a priori -> 3 failed to converge -> **912 fitted**.
**1477** interaction fields (not 1481; four duplicated input files were read alongside their
originals). **26,000** bins, a 200 x 130 grid at bin 50. Table 1 reproduces exactly from the
archived per-feature scores. Controls: spatial permutation test observed ELSA(d50) 0.075 with
20 domains against a 49-permutation null averaging 0.470 with 4-12 domains, p = 0.020;
concordance ARI 0.401, NMI 0.604, size-weighted purity 77.6%, p = 0.002.

Fuller detail is in project memory (`plosone-verified-counts`) and
`notation-to-code-groundwork.md`. **The pairing source matters**: 1477 comes from
`fantom-hierarchy.Robj`, not from `NICHES::ncomms8866_mouse`, which gives 1491 and will look
like a discrepancy if you use it by mistake.

## Open scientific threads, deliberately NOT part of this revision

All in `../../NEXT-PAPER-DIRECTIONS.md` at the repository root. The live one:

**Standardising before convolution recovers the ligand-receptor pairing signal; raw counts do
not.** Effective rank of the correlation matrix across the 1477 convolved fields, against an
exclusion-aware scrambled null: raw/minimum z = -0.51, raw/product z = -1.49, standardised/
minimum **z = +9.69**, standardised/product z = +2.20. This exposes a tension the paper does
not acknowledge, since it argues for raw counts and picks common minimum on cluster coherence.
Two cautions: the sign is opposite to what was predicted when the statistic was chosen and the
direction needs independent confirmation; and 8 scrambles floors the permutation p at 1/9, so
the z-score is the only informative summary.

Also recorded there: a naive receptor permutation recreates ~5% real FANTOM5 pairs, which
contaminated an earlier version of this test *and* the published Control A.

## Working notes

- R lives in the `emacs-r` distrobox, LaTeX in `bittensor-ubuntu`. Neither is on the host.
- Two long-running R jobs were killed during this session. **A tool timeout detaches the shell
  but does not kill R**; check `pgrep -af Rscript` before assuming a job died, or it will
  compete for CPU invisibly.
- Write R progress with `message()`, not `cat()`. `Rscript` block-buffers stdout to a file, so
  `cat()` output does not appear until the process exits and a running job looks frozen.
- The two `.DS_Store` deletions in `git status` predate this work. Leave them unstaged.
- Two untracked WIP files sit in `supp-4-lik-comparison/` (`04-run-diffusion.R`,
  `diffusion-function-testing.R`). They were deliberately not committed and still carry the
  `asp` bug fixed everywhere else.
