# Building the revision (PLOS ONE PONE-D-26-04349)

PLOS requires three uploads: a **Response to Reviewers**, a **Revised Manuscript with
Track Changes**, and a clean **Manuscript**. This file records how each is produced.

Everything builds inside the `bittensor-ubuntu` distrobox, which has TeX Live and
`latexdiff` (`/usr/bin/latexdiff`).

## Edit the manuscript as ordinary LaTeX

`plos-port/paper-plos.tex` carries **no track-changes markup**. Do not wrap edits in
`\add{}`, `\del{}`, or similar; the `oztrackchanges` package was removed on 2026-07-27.
The marked-up copy is generated mechanically by `latexdiff` at the end, which is complete
by construction and cannot miss an edit you forgot to mark.

## 1. Clean manuscript

```bash
distrobox enter bittensor-ubuntu -- bash -c \
  "cd '<repo>/Submission_PLOSOne/submission-2/plos-port' && \
   pdflatex -interaction=nonstopmode paper-plos.tex && \
   bibtex paper-plos && \
   pdflatex -interaction=nonstopmode paper-plos.tex && \
   pdflatex -interaction=nonstopmode paper-plos.tex"
```

Expected as of the baseline: 19 pages, 0 errors, 0 undefined references or citations.
`bibtex` reports one harmless "Repeated entry" for a duplicated `Granados2024` key in
`references-msbr.bib`; bibtex skips the second copy.

`paper-plos.pdf` is gitignored. Regenerate it rather than committing it.

## 2. Marked-up manuscript (track changes)

The baseline is the tag **`rev1-baseline`**, which points at the manuscript ported to the
PLOS class but with no revision edits applied.

```bash
cd <repo>/Submission_PLOSOne/submission-2/plos-port
git show rev1-baseline:Submission_PLOSOne/submission-2/plos-port/paper-plos.tex > /tmp/old.tex

distrobox enter bittensor-ubuntu -- bash -c \
  "cd '$PWD' && latexdiff --math-markup=off /tmp/old.tex paper-plos.tex > paper-plos-diff.tex && \
   pdflatex -interaction=nonstopmode paper-plos-diff.tex && \
   bibtex paper-plos-diff && \
   pdflatex -interaction=nonstopmode paper-plos-diff.tex && \
   pdflatex -interaction=nonstopmode paper-plos-diff.tex"
```

Deletions render struck through in red, insertions underlined in blue.

**Use `--math-markup=off`, not `--math-markup=whole`.** The R2.1 revision rewrote the
`alignedat` body of `eq:one-feat-model` (the process stage became additive, and a new
symbol `u_i(s)` was introduced). With `--math-markup=whole`, `latexdiff` still exits 0 but
emits markup inside that display which does not compile: the build produced **165 errors**
when checked on 2026-07-27. With `--math-markup=off` the same comparison compiles to 27
pages with **0 errors**. The cost is that changes *within* math displays are not
individually marked; the surrounding prose, which is where nearly all the revision lives,
is marked normally. Note the consequence for the reviewer-facing PDF: the revised model
equation will not be visibly flagged as changed, so the R2.1 response text should describe
the equation change in words, which it now does.

**Do not diff against the round-1 source in `submission-1/`.** That file uses the Elsevier
`elsarticle` class, so the diff would be dominated by the class port rather than by the
revision, and would be useless to a reviewer.

If a heavily rewritten table, equation, or figure environment produces output that will not
compile, the usual escapes are `--math-markup=off`, `--append-safecmd=<cmd>` for a macro
latexdiff is breaking, and `--exclude-textcmd=<cmd>` to stop it marking inside a command.
Reach for these only when a build actually fails.

Note that the colored boxes around citation numbers, equation numbers, and page numbers are
`hyperref` link decorations from the PLOS template. They appear in the clean build too and
are not diff markup.

## 2b. S9 Table (Supporting Information)

`plos-port/S9-notation-to-code.tex` is a standalone landscape document mapping each
mathematical object in the manuscript to the script and line that implements it (R2.2). It
has no bibliography and needs no `bibtex`:

```bash
distrobox enter bittensor-ubuntu -- bash -c \
  "cd '<repo>/Submission_PLOSOne/submission-2/plos-port' && \
   pdflatex -interaction=nonstopmode S9-notation-to-code.tex && \
   pdflatex -interaction=nonstopmode S9-notation-to-code.tex"
```

Expected: 4 pages, 0 errors. Upload the resulting PDF as the S9 Table Supporting
Information file. It is deliberately excluded from the `latexdiff` comparison, since it is
new in this revision and has no round-1 counterpart.

## 3. Response to Reviewers

```bash
distrobox enter bittensor-ubuntu -- bash -c \
  "cd '<repo>/Submission_PLOSOne/submission-2' && \
   pdflatex -interaction=nonstopmode response-to-reviewers.tex && \
   pdflatex -interaction=nonstopmode response-to-reviewers.tex"
```

The letter opens with an **"Open items" section listing every outstanding task**, including
the ones that have no `\todo` marker of their own, each cross-referenced to the relevant
response subsection and to the manuscript section, figure, or line it concerns. That is the
single place to look for what is left; start there.

### Status vocabulary: three words, used for both the item and the marker

"TODO" used to mean two different things: the inline red marker meant "something remains",
while `\status{TODO}` meant "not started", so an item could be `DONE` and still carry a red
`[TODO]`. There is now one vocabulary, three states, used for both:

| State | Meaning |
|---|---|
| `DONE` | finished and verified; nothing remains |
| `INCOMPLETE` | started, with real work already in the manuscript, but not finished |
| `TODO` | not started |

`\status{...}` labels a reviewer item as a whole. `\incomplete{...}` and `\todo{...}` mark one
specific outstanding task inside an item, and always carry the same word as that item's
status, so the two cannot disagree. A `DONE` item carries no marker.

Run

```bash
python3 check-response-letter.py
```

to verify: the closed vocabulary, the invariant, that the marker word matches the item status,
and that each subsection's status agrees with its row in the summary table. It exits non-zero
on any problem. All of those had drifted silently before the check existed.

**Before submitting**, set `\draftmodefalse` in the preamble. That removes the internal
banner, the Open-items section, and the status-summary table in one step (verified
2026-07-27: 21 pages drops to 17, with zero internal content in the output). Note that the
`\action` and `\authors` markers deliberately still render after the switch, as a loud
safety net, so also confirm none remain:

```bash
python3 check-response-letter.py    # also reports how many markers remain, and where
```

Do not try to count the markers with a bare `grep`: the obvious pattern silently matches
nothing (the leading backslash-a is consumed as an escape, so a broken gate looks like a
passing one), and a fixed-string grep additionally counts the two `\long\def` definitions in
the preamble. The checker splits the file by section, so it counts only real markers.


The `\todo` and `\incomplete` macros are defined `\long` so a marker may span paragraphs.
`\textcolor` and `\textbf` are not long, so a blank line inside a marker will still break the build with
"Paragraph ended before \@textcolor was complete."

## Figures

Updated 2026-10-08. One source, two builds, switched by `\embedfigstrue` / `\embedfigsfalse`
near the top of `paper-plos.tex`:

- **Review build (default, `\embedfigstrue`).** Figures are embedded so co-authors can read a
  single PDF. Figs 5 and 6 are scaled to share a page with their long captions; this affects
  the review PDF only.
- **Upload build (`\embedfigsfalse`).** Only the captions remain in the text, as PLOS requires,
  and the figures are uploaded separately as TIFFs.

Make the upload TIFFs, from `Submission_PLOSOne/submission-2/`:

```bash
bash export-submission-figures.sh      # writes submission-figures/Fig1.tif ... Fig6.tif
```

It fits each figure inside the PLOS limits (7.5 x 8.75 in, 300-600 dpi, 10 MB, checked against
journals.plos.org/plosone/s/figures on 2026-10-08) without upscaling. Figure sources:

| Fig | source | made by |
|-----|--------|---------|
| 1, 2 | `plos-port/img/` | original images, hand made |
| 3, 4 | `figures-regenerated/Fig3-composite.png`, `Fig4-composite.png` | `compose-figs-3-4.R` |
| 5 | `figures-regenerated/Fig5-composite-transposed.png` | `regenerate-fig5.R`, `fig5-input-panels.R`, `fig5-summary-panels.R`, then `compose-fig-5-transposed.py` |
| 6 | `figures-regenerated/Fig6-composite.png` | `regenerate-fig6.R` (reuses the 3D surfaces from `Fig6-relabeled.png`) |
| S10 | `figures-regenerated/S10_Fig_controls.png` | `make-controls-figure.R` |
| S11 | `figures-regenerated/S11_Fig_plugin-bias.png` | `code-spatial-smoothing/uncertainty/jensen-bias-figure.R` |

Figs 3 to 6 and S11 render text at 8 pt or larger at printed size, in Liberation Sans (metrically
identical to Arial). Figs 1 and 2 have not been checked for this.
