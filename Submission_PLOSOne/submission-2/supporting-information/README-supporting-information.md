# Supporting Information: what exists, what can be uploaded, what cannot

Prepared 2026-07-27 for reviewer requirement JR-5 ("Please upload a copy of Supporting
Information ... which you refer to in your text on page 19").

## The short version

**S1 through S8 cannot be uploaded to PLOS as Supporting Information files.** They total
**53.5 GB across 12,275 files**, and the largest single file is 11.6 GB. These are research
data archives, not supplementary documents, and they exceed any journal per-file limit by
three to four orders of magnitude. This is presumably why they were placed on FigShare in the
first place.

**S9 Table can and should be uploaded**, as a genuine 132 KB PDF. It is the only item here
that is a supplementary *document* rather than a data archive.

A decision is therefore needed from the authors, and it is a submission-strategy decision
rather than a technical one. The options are set out at the bottom of this file.

## Inventory, measured not estimated

| Item | Title | Size | Files | Largest file | Mostly |
|---|---|---|---|---|---|
| S1 File | Raw-to-bin | 3.6 GB | 9 | 3.7 GB | PDF (8), CSV (1) |
| S2 File | Feature refining | 1.3 GB | 2,926 | 569 MB | PNG (2,919) |
| S3 File | Raw vs norm, selected features | 37 MB | 24 | 17 MB | PNG (23) |
| S4 File | Model comparison, all features | 2.6 GB | 3,648 | 2.2 MB | PNG (2,736), CSV (912) |
| S5 File | HD embeddings | 38 GB | 1,547 | 11.6 GB | CSV (1,481), Robj (33) |
| S6 File | ELSA analysis | 322 MB | 1,330 | 96 MB | PNG (1,313), PDF (14) |
| S7 File | Clusters vs features | 2.2 GB | 1,310 | 1.3 GB | PNG (1,303) |
| S8 File | Morphogen fields | 5.5 GB | 1,481 | 3.8 MB | HTML (1,481) |
| **S9 Table** | **Notation-to-code mapping** | **132 KB** | **1** | **132 KB** | **PDF** |

Total S1-S8: 53.5 GB, 12,275 files.

## Current FigShare locations

These are the private-link URLs currently cited in the manuscript. Note that
`figshare.com/s/<hash>` links are *private shares*, not published DOIs. They are appropriate
for review but should be converted to public, DOI-minted items before publication, otherwise
the links can be revoked or expire and the citation becomes dead.

| Item | Link |
|---|---|
| S1 File | https://figshare.com/s/e921e8e36296cd8c57da |
| S2 File | https://figshare.com/s/6656c05d4552a59afa05 |
| S3 File | https://figshare.com/s/4885d58ba45b2f2e7689 |
| S4 File | https://figshare.com/s/9d266d5f0fd39632986d |
| S5 File | https://figshare.com/s/7bbe3feb729a9731acfa |
| S6 File | https://figshare.com/s/4867799e0e8b68fba5a8 |
| S7 File | https://figshare.com/s/2e305ce9ee20efd7ee47 |
| S8 File | https://figshare.com/s/2f66def6d22f1549b312 |

## What is in this directory

- `S9_Table.pdf` — the notation-to-code mapping, ready to upload as Supporting Information.
  Built from `../plos-port/S9-notation-to-code.tex`.

Nothing else is staged here, because nothing else is small enough to be worth staging.

## Where each item is cited in the manuscript

Confirmed against `../plos-port/paper-plos.tex` on 2026-07-27. S1 and S9 citations were added
during this revision; S3 through S8 were reassigned, because the original text cited
"Supplement 3/4/5" under a numbering that did not match the S1-S8 list.

| Item | Cited at |
|---|---|
| S1 File | Dataset preprocessing, binning paragraph |
| S2 File | Dataset preprocessing, feature-selection threshold |
| S3 File | Results, raw vs count-normalized; and the latent total-count paragraph |
| S4 File | Results, likelihood comparison, all-feature plots |
| S5 File | Results, 32 clustered spatial assays |
| S6 File | Results, convolution operator comparison, full results |
| S7 File | Results, feature-depth sweep, full results |
| S8 File | Results, interactive 3D interaction fields |
| S9 Table | Materials and methods, end of the implementation paragraph |

## The decision the authors need to make

1. **Upload S9 Table only, and keep S1-S8 on FigShare.** Least work. Requires converting the
   FigShare private links to public DOIs, and telling the editor plainly that the remaining
   items are 53 GB of derived data that belong in a repository. Most journals accept this.

2. **Additionally build small "representative" SI documents.** For example, S3 is already only
   37 MB and 23 PNGs, which could become a single modest PDF; S6's 14 summary PDFs are the
   interpretable part of its 1,330 files. This gives a reviewer something to open without a
   download, at the cost of curation work. Choosing what represents each archive is an author
   judgement, which is why it has not been done here.

3. **Move everything to a DOI-minting repository** (FigShare public, Zenodo, or Dryad) and
   cite DOIs throughout. Best for permanence and for the R2.9 reproducibility comment, and it
   is the option a reviewer checking reproducibility would prefer.

Whichever is chosen, the manuscript's Supporting information section and the response letter's
JR-5 both need to match it. They currently describe option 1.
