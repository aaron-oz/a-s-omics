# Review copies: start here

Snapshot PDFs of the current revision, committed so co-authors can open them directly from
GitHub without setting up LaTeX. Files 1-3 were rebuilt 2026-10-08 (file 4 is unchanged since 2026-07-30), with zero errors
and zero undefined references.

Read them in this order:

| File | Pages | What it is |
|---|---|---|
| `1_Revised-Manuscript_clean.pdf` | 32 | The revised manuscript as it would be read |
| `2_Revised-Manuscript_track-changes.pdf` | 34 | The same, versus the round-1 submission. Deletions struck through in red, insertions underlined in blue |
| `3_Response-to-Reviewers.pdf` | 21 | Point-by-point response. The "Open items" section near the start lists what is still outstanding |
| `4_S9-Table_notation-to-code.pdf` | 4 | New supplementary table mapping each mathematical symbol to the script and line implementing it (R2.2) |

If you are picking this work up rather than reviewing it, start instead with
`../REVISION-SUMMARY-2026-10-08.org`, a plain-language account of all the revision work and
where it has landed.

Also worth a look, outside this directory:

- `../figures-regenerated/` - Figures 3 to 6 and the graphical abstract, each built by a
  script there, at correct 1:1 aspect and with 8 pt or larger type at printed size.
- `../supporting-information/` - S9 Table, S10 Fig and S11 Fig, staged for upload.
- `../supporting-information/README-supporting-information.md` - measured inventory of S1
  through S9 and the Supporting Information decision that needs making.

## Two caveats when reading the track-changes copy

1. **Changes inside displayed equations are not marked.** The revision rewrote the model
   specification, but `latexdiff --math-markup=off` is required for the file to compile at all
   (the default `--math-markup=whole` produces 165 LaTeX errors on this source), and that
   switch suppresses markup within math displays. The equation change is described in words
   under R2.1 in the response letter.

2. **Figures 5 and 6 are scaled down in these PDFs** so each shares a page with its long
   caption. This affects the review PDFs only; the journal receives full-size TIFFs.

## Screen quality, NOT submission quality

These have been downsampled with ghostscript (`-dPDFSETTINGS=/ebook`) to bring the two
manuscript PDFs from 27 MB each to about 1 MB, so they are quick to open and do not bloat the
repository. Figures are legible on screen but the embedded images are no longer at full
resolution. **Do not upload these to the journal.** Build the submission copies from source per
`../README-build.md`, and remember that PLOS wants the figures as separate full-resolution
files anyway.

## These are snapshots, not the source of truth

The `.tex` files are. If you edit the sources, rebuild per `../README-build.md` and refresh
these copies, otherwise they will drift. They are committed for convenience of review, which
is why the build outputs in `../plos-port/` remain gitignored.
