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
  "cd '$PWD' && latexdiff --math-markup=whole /tmp/old.tex paper-plos.tex > paper-plos-diff.tex && \
   pdflatex -interaction=nonstopmode paper-plos-diff.tex && \
   bibtex paper-plos-diff && \
   pdflatex -interaction=nonstopmode paper-plos-diff.tex && \
   pdflatex -interaction=nonstopmode paper-plos-diff.tex"
```

Deletions render struck through in red, insertions underlined in blue. Verified working
on 2026-07-27 against representative edits (a scattered number correction and an inserted
clause): `latexdiff` exited 0 and the result compiled to 19 pages with no errors.

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

## 3. Response to Reviewers

```bash
distrobox enter bittensor-ubuntu -- bash -c \
  "cd '<repo>/Submission_PLOSOne/submission-2' && \
   pdflatex -interaction=nonstopmode response-to-reviewers.tex && \
   pdflatex -interaction=nonstopmode response-to-reviewers.tex"
```

**Before submitting**, set `\draftmodefalse` in the preamble to remove the internal banner
and the status-summary section, and confirm no `\todo` markers remain:

```bash
grep -c '\\todo{' response-to-reviewers.tex   # must be 0
```

The `\todo` macro is defined `\long` so a marker may span paragraphs. `\textcolor` and
`\textbf` are not long, so a blank line inside a `\todo{}` will still break the build with
"Paragraph ended before \@textcolor was complete."

## Figures

Figures are currently embedded in `paper-plos.tex` for a reviewable PDF. For the PLOS
upload they must be submitted as separate files and stripped from the `.tex`. They also
need regenerating for R2.6 (aspect-ratio distortion in Figures 3 and 4) and R1.3 (size and
organization), so do the stripping after that work, not before.
