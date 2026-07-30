#!/usr/bin/env bash
## Stream the 11 GB raw Stereo-seq GEM and keep only (a) a fixed panel of canonical
## anatomical marker genes and (b) reads falling inside the sub.bound FOV modeled in
## the paper. Produces controls/fov-markers-raw.tsv (~158 MB), the input used to build
## the finer-grained "MOSTA-plus" annotation (heart chambers, liver compartments).
##
## The GEM is a plain tab-separated table (geneID x y MIDCount) preceded by six
## '#'-comment header lines; coordinates are in the same raw frame as our grid, so no
## transformation is needed. awk is used rather than R because the full file does not
## fit comfortably in memory.
##
## usage: bash extract-fov-markers.sh
set -euo pipefail
cd "$(dirname "$0")/../.."

GEM="data-inputs/mouse-embryo-raw/full-raw-t16-data.tsv"
OUT="data-outputs/mouse-embryo-raw/controls/fov-markers-raw.tsv"

## sub.bound FOV, identical to supp-4 00-launch and control-B-extract-input.R
XMIN=15000; XMAX=25000; YMIN=11500; YMAX=18000

## Canonical markers, grouped by the structure they were chosen to resolve.
MARKERS="Myl2 Myl3 Myl4 Myl7 Sln Nppa Irx4 Hey2 Nr2f2 Kcnj3 Myh3 Myod1 Myog \
Alb Afp Hnf4a Apoa1 Apob Hba-a1 Hbb-bs Alas2 Gata1 Klf1 \
Shh Olig2 Nkx6-1 Nkx2-2 Foxa2 Mnx1 Isl1 Sim1 Pax3 Pax7 Msx1 Zic1 Wnt1 Lmx1a Gdf7 Lhx2 Lhx9 \
Nkx2-1 Sftpc Sox9 Acan Col2a1 Col1a1 Col3a1 Bmp4 Krt5 Krt14"

mkdir -p "$(dirname "$OUT")"
awk -v markers="$MARKERS" \
    -v xmin="$XMIN" -v xmax="$XMAX" -v ymin="$YMIN" -v ymax="$YMAX" '
BEGIN {
  FS = OFS = "\t"
  n = split(markers, m, /[ \t]+/)
  for (i = 1; i <= n; i++) keep[m[i]] = 1
  print "geneID", "x", "y", "MIDCount"
}
/^#/       { next }                      # GEM comment header
$1 == "geneID" { next }                  # column header
keep[$1] && $2 >= xmin && $2 <= xmax && $3 >= ymin && $3 <= ymax { print $1, $2, $3, $4 }
' "$GEM" > "$OUT"

echo "wrote $OUT ($(wc -l < "$OUT") lines, $(du -h "$OUT" | cut -f1))"
