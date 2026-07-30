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

### The two markup systems, which are deliberately different

The letter tracks two separate things, and they used to both be called "TODO", which was
ambiguous. They are now distinct:

- **`\status{...}`** is the state of a whole reviewer item, and is exactly one of
  `COMPLETE`, `PARTIAL`, or `NOT STARTED`. Nothing else is allowed.
- **`\action{...}`** and **`\authors{...}`** are individual outstanding tasks within an item.
  `\action` is a task anyone competent could do; `\authors` is a fact, confirmation, or
  decision only the authors can supply. They render as red `[ACTION: ...]` and
  `[AUTHORS: ...]`.

The invariant is that an item is `COMPLETE` if and only if it carries neither marker. Run

```bash
python3 check-response-letter.py
```

to verify all three properties at once: closed status vocabulary, the invariant, and
agreement between each subsection's status and its row in the summary table. It exits
non-zero on any problem. These three had all silently drifted before the check existed.

**Before submitting**, set `\draftmodefalse` in the preamble. That removes the internal
banner, the Open-items section, and the status-summary table in one step (verified
2026-07-27: 21 pages drops to 17, with zero internal content in the output). Note that the
`\action` and `\authors` markers deliberately still render after the switch, as a loud
safety net, so also confirm none remain:

```bash
# Use -F (fixed strings). Regex forms silently match nothing here, because the
# backslash-a is consumed as an escape, so a broken gate looks like a passing one.
# Verified 2026-07-27: returns 8 with markers present.
grep -c -F -e '\action{' -e '\authors{' response-to-reviewers.tex   # must be 0
```

The `\action` and `\authors` macros are defined `\long` so a marker may span paragraphs.
`\textcolor` and `\textbf` are not long, so a blank line inside a marker will still break the build with
"Paragraph ended before \@textcolor was complete."

## Figures

Figures are currently embedded in `paper-plos.tex` for a reviewable PDF. For the PLOS
upload they must be submitted as separate files and stripped from the `.tex`. They also
need regenerating for R2.6 (aspect-ratio distortion in Figures 3 and 4) and R1.3 (size and
organization), so do the stripping after that work, not before.

**Figures 5 and 6 are deliberately scaled down in the embedded PDF, and this is temporary.**
Both source images are tall portraits (5100 x 6600), and the R2.7 rewrite made both captions
much longer. Together they exceeded `\textheight`, so LaTeX floated them past the
bibliography to pages 26-27, ten pages after their citations, with the caption text running
off the bottom of the page. The fix applied is `[p]` placement, `\captionsetup{font=small}`
on those two captions, and `width=0.705`/`0.745\linewidth`; they now sit on pages 16 and 17,
immediately after the text that cites them, with complete captions.

The scaling works against R1.3, which asks for *larger* figures, so do not treat it as the
answer to that comment. It is a stopgap for the reviewable PDF only. Two things dissolve it:
the figures are being regenerated and reorganized for R1.3 (a less tall aspect ratio would
remove the constraint entirely), and in the actual PLOS submission the figures are uploaded
as separate full-resolution files with captions set in the text, so the float-fitting problem
does not arise. When you strip the figures for upload, drop the `[p]`, the `\captionsetup`,
and the width scaling along with them.
