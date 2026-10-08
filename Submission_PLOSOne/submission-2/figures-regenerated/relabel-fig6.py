#!/usr/bin/env python3
"""Correct the baked-in "1481" in Fig 6's full-depth row to "1477".

Fig 6 has no generating layout script (it was assembled by hand), so this edits the raster
directly. The full-depth row carries "1481" in four places: the UMAP panel title, the spatial
panel title, and both UMAP axis labels ("UMAP1481feature_1/_2"). Only the digits "81" are
repainted as "77"; Helvetica digits share one advance width, so nothing else moves.

Why 1477 and not 1449 (checked 2026-10-08 against m-est-min-matrix.rds, which is the
1481-row clustered object): the run clustered 1481 columns = all 1477 FANTOM5 mechanisms
plus four duplicated input columns (PDGFC/PDGFD x PDGFRA/PDGFRB). The other rows' labels
("10 Features" ... "1300 Features") count mechanisms entered, so the full-depth row is
"all 1477 mechanisms". 1449 is a different quantity (numerically distinct fields: the
common minimum maps 28 further mechanisms onto fields identical to another mechanism's),
and belongs in the caption, not the label.

Usage (no system Python packages needed):
  uv run --no-project --with pillow --with numpy python relabel-fig6.py IN.png OUT.png
IN must be the original, unrelabeled image, e.g.
  git show "rev1-baseline:Submission_PLOSOne/submission-2/plos-port/img/Spatial Domains.png" > IN.png
Output was copied over plos-port/img/Spatial Domains.png on 2026-10-08.
"""
import sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont

FONT = "/usr/share/fonts/liberation-sans-fonts/LiberationSans-Regular.ttf"

# (x0, y0, x1, y1) boxes holding exactly one text line, the rotation (degrees CCW that
# makes the text read left to right), and which glyphs (0-based) are the "8" and final "1".
LABELS = [
    ((95, 3694, 420, 3730), 0, (2, 3)),     # "1481 Features | UMAP | Resolution = 1"
    ((1495, 3700, 1840, 3746), 0, (2, 3)),  # "1481 Features | Spatial | Resolution = 1"
    ((540, 4548, 820, 4580), 0, (6, 7)),    # "UMAP1481feature_1"
    ((14, 3950, 46, 4500), 90, (6, 7)),     # "UMAP1481feature_2", rotated
]


def glyph_spans(mask):
    cols = mask.any(axis=0)
    spans, start = [], None
    for i, c in enumerate(cols):
        if c and start is None:
            start = i
        elif not c and start is not None:
            spans.append((start, i - 1)); start = None
    if start is not None:
        spans.append((start, len(cols) - 1))
    return spans


def fix(img, box, rot, idx):
    x0, y0, x1, y1 = box
    crop = img.crop(box)
    if rot:
        crop = crop.rotate(-rot, expand=True)  # make it read left to right
    a = np.asarray(crop.convert("L"))
    mask = a > 110
    spans = glyph_spans(mask)
    s8, s1 = spans[idx[0]], spans[idx[1]]
    one = spans[idx[0] - 2]  # the leading "1", for the digit height
    rows = np.where(mask[:, one[0]:one[1] + 1].any(axis=1))[0]
    top, bot = rows.min(), rows.max()
    h = bot - top + 1
    draw = ImageDraw.Draw(crop)
    draw.rectangle([s8[0] - 1, top - 2, s1[1] + 1, bot + 2], fill=(0, 0, 0))
    # size the font so a digit's ink height matches the original digits
    size = h
    for _ in range(20):
        f = ImageFont.truetype(FONT, size)
        bb = f.getbbox("7")
        if bb[3] - bb[1] >= h:
            break
        size += 1
    f = ImageFont.truetype(FONT, size)
    # advance per digit, measured from the original "4" to "8"
    adv = s8[0] - spans[idx[0] - 1][0]
    bb7 = f.getbbox("7")
    for k in range(2):
        x = s8[0] + k * adv - bb7[0]
        draw.text((x, bot - bb7[3] + 1), "7", font=f, fill=(255, 255, 255))
    if rot:
        crop = crop.rotate(rot, expand=True)
    img.paste(crop, (x0, y0))


def main(src, dst):
    img = Image.open(src).convert("RGB")
    for box, rot, idx in LABELS:
        fix(img, box, rot, idx)
    img.save(dst, dpi=img.info.get("dpi", (600, 600)))


if __name__ == "__main__":
    main(*sys.argv[1:3])
