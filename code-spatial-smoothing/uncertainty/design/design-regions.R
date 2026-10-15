## Frozen regions and slices for the uncertainty pivot (Phase 2, PIVOT-PLAN-2026-10-14.org).
## Single source of truth: every later script reads design-regions.rds, never re-derives these.
## Chosen from anatomy (MOSTA annotation) and total counts only; no gene, signaling field or
## posterior quantity was looked at when choosing them (design-maps.R shows what was looked at).
## Bin index: i = (x - 15025)/50 + 1 in 1..200, j = (y - 11525)/50 + 1 in 1..130.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript design-regions.R
suppressMessages({library(data.table); library(ggplot2)})
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
a <- fread(file.path(repo, "data-outputs/mouse-embryo-raw/controls/e16-fov-annotation.csv"))
a[, `:=`(i = as.integer((x - 15025)/50 + 1), j = as.integer((y - 11525)/50 + 1))]
setkey(a, i, j)

## ---- fixed boxes: 15 x 15 bins (225 bins), roughly 1.5 fitted correlation lengths on a side
boxes <- data.table(
  name  = c("box-heart-1", "box-heart-2", "box-liver", "box-cord-1", "box-cord-2"),
  i0 = c(100, 115, 170,  30,  85), j0 = c( 85, 105, 100,   6,   3))
boxes[, `:=`(i1 = i0 + 14L, j1 = j0 + 14L)]
## ---- anatomical regions: every bin carrying the MOSTA label
anat <- c("Heart", "Liver", "Lung", "Spinal cord", "Thymus")

regions <- list()
for (r in seq_len(nrow(boxes))) {
  b <- boxes[r]
  regions[[b$name]] <- a[i >= b$i0 & i <= b$i1 & j >= b$j0 & j <= b$j1, .(x, y, i, j, annotation)]
}
for (nm in anat) regions[[paste0("anat-", tolower(gsub(" ", "-", nm)))]] <- a[annotation == nm, .(x, y, i, j, annotation)]

## ---- slices: straight transects, one bin per unit step along the longer axis
slices <- data.table(
  name = c("slice-cord-to-muscle", "slice-heart-to-liver", "slice-along-cord", "slice-thymus-to-heart"),
  ia = c( 50,  85,  10,  40), ja = c( 1, 95, 10, 75),
  ib = c( 50, 200, 150, 140), jb = c(51, 95,  5, 75))
slice_bins <- list()
for (r in seq_len(nrow(slices))) {
  s <- slices[r]; n <- max(abs(s$ib - s$ia), abs(s$jb - s$ja)) + 1
  t <- seq(0, 1, length.out = n)
  ii <- as.integer(round(s$ia + t * (s$ib - s$ia))); jj <- as.integer(round(s$ja + t * (s$jb - s$ja)))
  d <- unique(data.table(pos = seq_len(n), i = ii, j = jj), by = c("i", "j"))
  slice_bins[[s$name]] <- a[d, on = .(i, j)][order(pos), .(pos, x, y, i, j, annotation)]
}

## ---- report composition
cat("Regions (bins, label composition):\n")
for (nm in names(regions)) {
  r <- regions[[nm]]; tab <- sort(table(r$annotation), decreasing = TRUE)
  cat(sprintf("  %-18s %5d bins  %s\n", nm, nrow(r),
              paste(sprintf("%s %.0f%%", names(tab)[1:min(3, length(tab))],
                            100 * tab[1:min(3, length(tab))] / nrow(r)), collapse = ", ")))
}
cat("Slices (bins, labels crossed in order):\n")
for (nm in names(slice_bins)) {
  s <- slice_bins[[nm]]; rl <- rle(s$annotation)
  cat(sprintf("  %-22s %3d bins: %s\n", nm, nrow(s),
              paste(sprintf("%s(%d)", rl$values, rl$lengths), collapse = " > ")))
}
saveRDS(list(boxes = boxes, anat = anat, regions = regions, slices = slices, slice_bins = slice_bins),
        "design-regions.rds")

## ---- overlay figure
big <- a[, .N, by = annotation][N >= 500, annotation]
a[, ann := ifelse(annotation %in% big, annotation, "other")]
sl <- rbindlist(slice_bins, idcol = "slice")
g <- ggplot(a, aes(i, j)) + geom_raster(aes(fill = ann), alpha = 0.6) + coord_fixed() +
  geom_rect(data = boxes, aes(xmin = i0 - 0.5, xmax = i1 + 0.5, ymin = j0 - 0.5, ymax = j1 + 0.5),
            inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.7) +
  geom_text(data = boxes, aes(x = (i0 + i1)/2, y = j1 + 3, label = name), size = 3) +
  geom_path(data = sl, aes(i, j, group = slice), color = "black", linewidth = 0.9) +
  geom_text(data = slices, aes(x = ia, y = ja + 2, label = name), size = 3, hjust = 0) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme_bw(9) + labs(fill = NULL, title = "Frozen boxes (rectangles) and slices (lines)")
ggsave("design-regions.png", g, width = 11, height = 7, dpi = 130)
