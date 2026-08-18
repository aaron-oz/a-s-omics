#!/usr/bin/env python3
"""Assemble Figure 5 (convolution operator comparison) from the regenerated panels.

Replaces hand assembly in Illustrator. Figure 5 is a portrait stack rather than a grid, so it
gets its own script: three operator columns, a top block of four example mechanisms each
showing ligand and receptor above their interaction field, then the resulting domains and the
matching entrogram, then two summary panels spanning the width.

The panels themselves come from regenerate-fig5.R, which applies the coord_fixed() correction
that Figure 5 needed (R2.6 as extended to this figure). Nothing is re-clustered anywhere.

Usage: python3 compose-fig-5.py [width_in] [dpi]
"""
import sys
from PIL import Image, ImageDraw, ImageFont

W_IN = float(sys.argv[1]) if len(sys.argv) > 1 else 7.5
DPI  = int(sys.argv[2]) if len(sys.argv) > 2 else 400
W    = int(W_IN * DPI)

P = "panels-fig5"
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

def load(name):
    return Image.open(f"{P}/{name}.png").convert("RGB")

def fit(im, w):
    return im.resize((w, round(im.height * w / im.width)), Image.LANCZOS)

# ---- geometry ------------------------------------------------------------------------
PAD      = round(0.030 * DPI)
GUT      = round(0.055 * DPI)          # gap between operator columns
LEFT     = round(0.30 * DPI)           # row-label gutter
COLW     = (W - 2 * PAD - LEFT - 2 * GUT) // 3
LR_W     = round(COLW * 0.30)          # small ligand / receptor panels
INT_W    = COLW - LR_W - round(0.02 * DPI)
F_HEAD   = font(round(0.105 * DPI), True)
F_SUB    = font(round(0.072 * DPI))
F_ROW    = font(round(0.078 * DPI), True)
F_MECH   = font(round(0.066 * DPI), True)
F_TITLE  = font(round(0.135 * DPI), True)

lr_h  = fit(load("Fig5_min_WNT4-FZD6_ligand"), LR_W).height
int_h = fit(load("Fig5_min_WNT4-FZD6_interaction"), INT_W).height
mech_h = max(2 * lr_h + round(0.012 * DPI), int_h) + round(0.075 * DPI)
map_h  = fit(load("Fig5_min_domains"), COLW).height
summ_h = fit(load("Fig5_summary_elsa-vs-paramD"), (W - 2 * PAD - GUT) // 2).height

TITLE_H = round(0.30 * DPI)
HEAD_H  = round(0.20 * DPI)
H = (TITLE_H + HEAD_H + 4 * mech_h + round(0.10 * DPI) + map_h + round(0.10 * DPI)
     + map_h + round(0.14 * DPI) + summ_h + 2 * PAD)

canvas = Image.new("RGB", (W, H), "white")
dr = ImageDraw.Draw(canvas)

def ctext(x, y, s, f, fill="black"):
    b = dr.textbbox((0, 0), s, font=f)
    dr.text((x - (b[2] - b[0]) / 2, y), s, font=f, fill=fill)

def rtext(x, y, s, f, fill="black"):
    b = dr.textbbox((0, 0), s, font=f)
    dr.text((x - (b[2] - b[0]), y), s, font=f, fill=fill)

def colx(i):
    return PAD + LEFT + i * (COLW + GUT)

# ---- title and operator headers ------------------------------------------------------
y = PAD
ctext(W / 2, y, "Comparison of ligand-receptor convolution operators", F_TITLE)
y += TITLE_H
for i, (_, label) in enumerate(OPS):
    ctext(colx(i) + COLW / 2, y, label, F_HEAD)
y += HEAD_H

# ---- top block: four mechanisms, ligand and receptor beside their interaction field ---
for m in MECHS:
    for i, (op, _) in enumerate(OPS):
        x = colx(i)
        canvas.paste(fit(load(f"Fig5_{op}_{m}_ligand"), LR_W), (x, y))
        canvas.paste(fit(load(f"Fig5_{op}_{m}_receptor"), LR_W), (x, y + lr_h + round(0.012 * DPI)))
        canvas.paste(fit(load(f"Fig5_{op}_{m}_interaction"), INT_W),
                     (x + LR_W + round(0.02 * DPI), y))
    dr.text((PAD, y + round(0.02 * DPI)), m, font=F_MECH, fill="black")
    y += mech_h

# ---- resulting domains, then the matching entrogram ----------------------------------
y += round(0.10 * DPI)
for i, (op, _) in enumerate(OPS):
    canvas.paste(fit(load(f"Fig5_{op}_domains"), COLW), (colx(i), y))
rtext(PAD + LEFT - round(0.04 * DPI), y + map_h / 2 - round(0.04 * DPI), "Domains", F_ROW)
y += map_h + round(0.10 * DPI)
for i, (op, _) in enumerate(OPS):
    canvas.paste(fit(load(f"Fig5_{op}_entrogram"), COLW), (colx(i), y))
rtext(PAD + LEFT - round(0.04 * DPI), y + map_h / 2 - round(0.04 * DPI), "Entrogram", F_ROW)
y += map_h + round(0.14 * DPI)

# ---- two summary panels spanning the width -------------------------------------------
sw = (W - 2 * PAD - GUT) // 2
canvas.paste(fit(load("Fig5_summary_elsa-vs-paramD"), sw), (PAD, y))
canvas.paste(fit(load("Fig5_summary_elsa-vs-numclust"), sw), (PAD + sw + GUT, y))

# ---- panel letters -------------------------------------------------------------------
letters = [("A", PAD + LEFT, PAD + TITLE_H + HEAD_H),
           ("B", PAD + LEFT, PAD + TITLE_H + HEAD_H + 4 * mech_h + round(0.10 * DPI)),
           ("C", PAD + LEFT, PAD + TITLE_H + HEAD_H + 4 * mech_h + round(0.20 * DPI) + map_h),
           ("D", PAD, y), ("E", PAD + sw + GUT, y)]
for s, lx, ly in letters:
    dr.text((lx - round(0.10 * DPI), ly), s, font=F_HEAD, fill="black")

# PLOS ONE caps a figure at 7.5 x 8.75 in. Figure 5 carries twenty panels in a portrait
# stack, so at full width it overruns the height limit. That is the same constraint that
# forced the published version to be scaled to 0.705\linewidth in the manuscript. Scale the
# finished composite to fit the box rather than silently exceeding it.
MAX_H_IN = 8.75
h_in = H / DPI
if h_in > MAX_H_IN:
    k = MAX_H_IN / h_in
    canvas = canvas.resize((round(W * k), round(H * k)), Image.LANCZOS)
    print(f"scaled by {k:.3f} to meet the {MAX_H_IN} in height limit")

canvas.save("Fig5-composite.png", dpi=(DPI, DPI))
print(f"Fig5-composite.png  {canvas.width / DPI:.2f} x {canvas.height / DPI:.2f} in @ {DPI} dpi"
      f"  ({canvas.width} x {canvas.height} px)")
