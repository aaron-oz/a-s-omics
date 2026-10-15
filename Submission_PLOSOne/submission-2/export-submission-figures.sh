#!/usr/bin/env bash
# Export the six main figures as PLOS-compliant TIFFs into submission-figures/.
# PLOS One (checked 2026-10-08, journals.plos.org/plosone/s/figures): TIFF or EPS; at most
# 7.5 x 8.75 in; 300-600 dpi; at most 10 MB; text 8-12 pt in Arial, Times or Symbol.
# Each figure is fitted inside 7.5 x 8.75 in at 600 dpi, then dropped to 400 dpi if the
# LZW-compressed file would exceed 10 MB.
# Run from Submission_PLOSOne/submission-2/:  bash export-submission-figures.sh
set -euo pipefail
out=submission-figures; mkdir -p "$out"
declare -A SRC=(
  [Fig1]="plos-port/img/Project Schematic.png"
  [Fig2]="plos-port/img/Data Parameters.png"
  [Fig3]="figures-regenerated/Fig3-composite.png"
  [Fig4]="figures-regenerated/Fig4-composite.png"
  [Fig5]="figures-regenerated/Fig5-composite-transposed.png"
  [Fig6]="figures-regenerated/Fig6-composite.png"
)
for f in Fig1 Fig2 Fig3 Fig4 Fig5 Fig6; do
  # never upscale: an image drawn at a lower resolution gets that resolution, so it still
  # prints at its intended size (Figs 3-4 are drawn 7.5 in wide at 400 dpi)
  read pw ph < <(magick identify -format "%w %h\n" "${SRC[$f]}")
  native=$(( pw * 10 / 75 )); hn=$(( ph * 100 / 875 )); [ "$hn" -gt "$native" ] && native=$hn   # the larger bound keeps both sides within the limit
  for dpi in 600 400 300; do
    [ "$dpi" -gt "$native" ] && dpi=$native
    w=$((75 * dpi / 10)); h=$((875 * dpi / 100))
    magick "${SRC[$f]}" -alpha remove -alpha off -resize "${w}x${h}>" \
      -units PixelsPerInch -density "$dpi" -compress lzw "$out/$f.tif"
    sz=$(stat -c %s "$out/$f.tif")
    if [ "$sz" -le 10000000 ]; then break; fi
  done
  identify -format "$f.tif  %wx%h px  %[fx:w/resolution.x]x%[fx:h/resolution.y] in  %x dpi  " "$out/$f.tif"
  echo "$((sz / 1000)) kB"
done
