#!/usr/bin/env python3
"""Assemble Figure 5 (convolution operator comparison) from the regenerated panels.

Replaces hand assembly in Illustrator. Figure 5 is a portrait stack rather than a grid, so it
gets its own script rather than sharing compose-figs-3-4.R.

The ligand and receptor fields are l.est and r.est, which do not depend on the convolution
operator. The published layout repeated them inside each of the three operator columns, so
eight panels were drawn three times each; verified byte-identical across the columns. They are
shown once here, in a left gutter beside the mechanism they belong to, and the width that
recovers goes to the interaction fields, which is what Reviewer 1 asked for in R1.3.

The panels themselves come from regenerate-fig5.R, which applies the coord_fixed() correction
Figure 5 needed (R2.6 as extended to this figure), and fig5-summary-panels.R. Nothing is
re-clustered.

Usage: python3 compose-fig-5.py [width_in] [dpi]
"""
import sys
from PIL import Image, ImageDraw, ImageFont

W_IN = float(sys.argv[1]) if len(sys.argv) > 1 else 7.5
DPI  = int(sys.argv[2]) if len(sys.argv) > 2 else 400
MAX_H_IN = 8.75                      # PLOS ONE figure box is 7.5 x 8.75 in
W    = int(W_IN * DPI)

P     = "panels-fig5"
OPS   = [("prod", "Product"), ("gm", "Geometric mean"), ("min", "Common minimum")]
MECHS = ["WNT4-FZD6", "RSPO1-LGR6", "DHH-PTCH1", "VEGFA-KDR"]

def font(px, bold=False):
    for path in ("/usr/share/fonts/dejavu-sans-fonts/DejaVuSans%s.ttf",
                 "/usr/share/fonts/truetype/dejavu/DejaVuSans%s.ttf"):
        try:
            return ImageFont.truetype(path % ("-Bold" if bold else ""), px)
        except OSError:
            continue
    return ImageFont.load_default()

def load(n): return Image.open(f"{P}/{n}.png").convert("RGB")
def fit(im, w): return im.resize((w, round(im.height * w / im.width)), Image.LANCZOS)

I = lambda x: round(x * DPI)         # inches to pixels

PAD, GUT = I(0.030), I(0.038)
LBL_W, LR_W, LR_GAP = I(0.30), I(0.52), I(0.030)
LEFT = LBL_W + LR_W + I(0.04)
COLW = (W - 2 * PAD - LEFT - 2 * GUT) // 3

F_TITLE, F_HEAD = font(I(0.135), True), font(I(0.100), True)
F_ROW, F_MECH   = font(I(0.075), True), font(I(0.068), True)

lr_h   = fit(load(f"Fig5_input_{MECHS[0]}_ligand"), LR_W).height
int_h  = fit(load(f"Fig5_min_{MECHS[0]}_interaction"), COLW).height
mech_h = max(int_h, 2 * lr_h + LR_GAP) + I(0.055)
map_h  = fit(load("Fig5_min_domains"), COLW).height
sw     = (W - 2 * PAD - GUT) // 2
summ_h = fit(load("Fig5_summary_elsa-vs-paramD"), sw).height

TITLE_H, HEAD_H = I(0.30), I(0.20)
H = (PAD + TITLE_H + HEAD_H + 4 * mech_h + I(0.10) + map_h + I(0.10) + map_h
     + I(0.14) + summ_h + PAD)

canvas = Image.new("RGB", (W, H), "white")
dr = ImageDraw.Draw(canvas)

def ctext(x, y, s, f, fill="black"):
    b = dr.textbbox((0, 0), s, font=f); dr.text((x - (b[2] - b[0]) / 2, y), s, font=f, fill=fill)

def rtext(x, y, s, f, fill="black"):
    b = dr.textbbox((0, 0), s, font=f); dr.text((x - (b[2] - b[0]), y), s, font=f, fill=fill)

colx = lambda i: PAD + LEFT + i * (COLW + GUT)

y = PAD
ctext(W / 2, y, "Comparison of ligand-receptor convolution operators", F_TITLE)
y += TITLE_H
ctext(PAD + LBL_W + LR_W / 2, y, "Inputs", F_HEAD, fill=(80, 80, 80))
for i, (_, label) in enumerate(OPS):
    ctext(colx(i) + COLW / 2, y, label, F_HEAD)
y += HEAD_H

top_y = y
for m in MECHS:
    lx = PAD + LBL_W
    canvas.paste(fit(load(f"Fig5_input_{m}_ligand"), LR_W), (lx, y))
    canvas.paste(fit(load(f"Fig5_input_{m}_receptor"), LR_W), (lx, y + lr_h + LR_GAP))
    for i, (op, _) in enumerate(OPS):
        canvas.paste(fit(load(f"Fig5_{op}_{m}_interaction"), COLW), (colx(i), y))
    rtext(PAD + LBL_W - I(0.045), y + I(0.02), m, F_MECH)
    y += mech_h

y += I(0.10)
dom_y = y
for i, (op, _) in enumerate(OPS):
    canvas.paste(fit(load(f"Fig5_{op}_domains"), COLW), (colx(i), y))
rtext(PAD + LEFT - I(0.045), y + map_h / 2 - I(0.04), "Domains", F_ROW)
y += map_h + I(0.10)
ent_y = y
for i, (op, _) in enumerate(OPS):
    canvas.paste(fit(load(f"Fig5_{op}_entrogram"), COLW), (colx(i), y))
rtext(PAD + LEFT - I(0.045), y + map_h / 2 - I(0.04), "Entrogram", F_ROW)
y += map_h + I(0.14)

canvas.paste(fit(load("Fig5_summary_elsa-vs-paramD"), sw), (PAD, y))
canvas.paste(fit(load("Fig5_summary_elsa-vs-numclust"), sw), (PAD + sw + GUT, y))

for s, lx, ly in [("A", PAD, top_y), ("B", PAD, dom_y), ("C", PAD, ent_y),
                  ("D", PAD, y), ("E", PAD + sw + GUT - I(0.10), y)]:
    dr.text((lx, ly), s, font=F_HEAD, fill="black")

h_in = H / DPI
if h_in > MAX_H_IN:
    k = MAX_H_IN / h_in
    canvas = canvas.resize((round(W * k), round(H * k)), Image.LANCZOS)
    print(f"scaled by {k:.3f} to meet the {MAX_H_IN} in height limit")

canvas.save("Fig5-composite.png", dpi=(DPI, DPI))
print(f"Fig5-composite.png  {canvas.width / DPI:.2f} x {canvas.height / DPI:.2f} in @ {DPI} dpi")
print(f"interaction panel width: {COLW * (canvas.width / W) / DPI:.2f} in")
