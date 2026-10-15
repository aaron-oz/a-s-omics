#!/usr/bin/env python3
"""Figure 5, transposed: operators as rows, mechanisms as columns.

Alternative to compose-fig-5.py, which keeps the published portrait arrangement. Same panels,
same data, nothing re-clustered. The point of comparison is page area: the portrait version is
height-bound and ends up 5.66 x 8.75 in, leaving 1.84 in of permitted width unused, while this
one runs full width.

Layout: a shared inputs row (each mechanism's ligand and receptor side by side), then one row
per convolution operator across the four mechanisms, then the resulting domains and entrograms
as their own rows, then the two summary panels.

Usage: python3 compose-fig-5-transposed.py [width_in] [dpi]
"""
import sys
from PIL import Image, ImageDraw, ImageFont

W_IN = float(sys.argv[1]) if len(sys.argv) > 1 else 6.69
DPI  = int(sys.argv[2]) if len(sys.argv) > 2 else 600
MAX_H_IN = 8.75
W = int(W_IN * DPI)

P     = "panels-fig5"
OPS   = [("prod", "Product"), ("gm", "Geometric\nmean"), ("min", "Common\nminimum")]
MECHS = ["WNT4-FZD6", "RSPO1-LGR6", "DHH-PTCH1", "VEGFA-KDR"]

def font(px, bold=False):
    # Liberation Sans is metrically identical to Arial, one of the PLOS figure fonts.
    # Fail loudly rather than fall back to PIL's tiny bitmap font.
    path = "/usr/share/fonts/liberation-sans-fonts/LiberationSans-%s.ttf" % ("Bold" if bold else "Regular")
    return ImageFont.truetype(path, px)

def load(n): return Image.open(f"{P}/{n}.png").convert("RGB")
def fit(im, w): return im.resize((w, round(im.height * w / im.width)), Image.LANCZOS)
I = lambda x: round(x * DPI)

PAD, GUT = I(0.030), I(0.032)
LEFT = I(0.62)                                     # operator row labels
COLW = (W - 2 * PAD - LEFT - 3 * GUT) // 4         # four mechanism columns
MAPW = (W - 2 * PAD - LEFT - 2 * GUT) // 3         # three operator columns for the map rows

PT8 = I(8 / 72)                                    # PLOS: text 8-12 pt in print
F_TITLE, F_HEAD = font(PT8, True), font(PT8, True)
F_ROW, F_SMALL  = font(PT8, True), font(PT8)

in_w   = (COLW - I(0.02)) // 2
in_h   = fit(load(f"Fig5_input_{MECHS[0]}_ligand"), in_w).height
int_h  = fit(load(f"Fig5_min_{MECHS[0]}_interaction"), COLW).height
map_h  = fit(load("Fig5_min_domains"), MAPW).height
sw     = (W - 2 * PAD - GUT) // 2
summ_h = fit(load("Fig5_summary_elsa-vs-paramD"), sw).height

TITLE_H, HEAD_H = 0, I(0.17)
H = (PAD + TITLE_H + HEAD_H + in_h + I(0.10) + 3 * (int_h + I(0.045)) + I(0.12)
     + 2 * (map_h + I(0.12)) + I(0.04) + summ_h + PAD)

canvas = Image.new("RGB", (W, H), "white")
dr = ImageDraw.Draw(canvas)

def ctext(x, y, s, f, fill="black"):
    b = dr.textbbox((0, 0), s, font=f); dr.text((x - (b[2] - b[0]) / 2, y), s, font=f, fill=fill)
def rtext(x, y, s, f, fill="black"):
    lines = s.split("\n"); lh = I(0.125)
    y0 = y - (len(lines) - 1) * lh / 2
    for k, ln in enumerate(lines):
        b = dr.textbbox((0, 0), ln, font=f); dr.text((x - (b[2] - b[0]), y0 + k * lh), ln, font=f, fill=fill)

colx  = lambda i: PAD + LEFT + i * (COLW + GUT)
mapx  = lambda i: PAD + LEFT + i * (MAPW + GUT)

y = PAD
y += TITLE_H
for i, m in enumerate(MECHS):
    ctext(colx(i) + COLW / 2, y, m, F_HEAD)
y += HEAD_H

# shared inputs: each mechanism's ligand and receptor, side by side, shown once
in_y = y
for i, m in enumerate(MECHS):
    canvas.paste(fit(load(f"Fig5_input_{m}_ligand"), in_w), (colx(i), y))
    canvas.paste(fit(load(f"Fig5_input_{m}_receptor"), in_w), (colx(i) + in_w + I(0.02), y))
rtext(PAD + LEFT - I(0.05), y + in_h / 2 - I(0.04), "Inputs", F_ROW, fill=(80, 80, 80))
y += in_h + I(0.10)

# one row per operator
op_y = y
for op, label in OPS:
    for i, m in enumerate(MECHS):
        canvas.paste(fit(load(f"Fig5_{op}_{m}_interaction"), COLW), (colx(i), y))
    rtext(PAD + LEFT - I(0.05), y + int_h / 2 - I(0.04), label, F_ROW)
    y += int_h + I(0.045)
y += I(0.12) - I(0.045)

# domains, then entrograms, one column per operator
dom_y = y
for i, (op, label) in enumerate(OPS):
    canvas.paste(fit(load(f"Fig5_{op}_domains"), MAPW), (mapx(i), y))
    ctext(mapx(i) + MAPW / 2, y - I(0.13), label.replace("\n", " "), F_SMALL, fill=(80, 80, 80))
rtext(PAD + LEFT - I(0.05), y + map_h / 2 - I(0.04), "Domains", F_ROW)
y += map_h + I(0.12)
ent_y = y
for i, (op, _) in enumerate(OPS):
    canvas.paste(fit(load(f"Fig5_{op}_entrogram"), MAPW), (mapx(i), y))
rtext(PAD + LEFT - I(0.05), y + map_h / 2 - I(0.04), "Entrogram", F_ROW)
y += map_h + I(0.12) + I(0.04)

canvas.paste(fit(load("Fig5_summary_elsa-vs-paramD"), sw), (PAD, y))
canvas.paste(fit(load("Fig5_summary_elsa-vs-numclust"), sw), (PAD + sw + GUT, y))

for s, lx, ly in [("A", PAD, in_y), ("B", PAD, op_y), ("C", PAD, dom_y),
                  ("D", PAD, ent_y), ("E", PAD + I(0.04), y + I(0.03)), ("F", PAD + sw + GUT + I(0.04), y + I(0.03))]:
    # E and F sit on the dark summary panels, so they are drawn in white
    dr.text((lx, ly), s, font=F_HEAD, fill="white" if s in "EF" else "black")

h_in = H / DPI
if h_in > MAX_H_IN:
    k = MAX_H_IN / h_in
    canvas = canvas.resize((round(W * k), round(H * k)), Image.LANCZOS)
    print(f"scaled by {k:.3f} to meet the {MAX_H_IN} in height limit")
canvas.save("Fig5-composite-transposed.png", dpi=(DPI, DPI))
fw, fh = canvas.width / DPI, canvas.height / DPI
print(f"Fig5-composite-transposed.png  {fw:.2f} x {fh:.2f} in @ {DPI} dpi")
print(f"interaction panel width: {COLW * (canvas.width / W) / DPI:.2f} in;  page area {fw*fh:.1f} sq in")
