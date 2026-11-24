## Blind version of threshold-examples.png for a second rater (REVIEW-2026-10-29.org item P4):
## the same five panels, same seed-drawn instances (threshold-examples.csv), same color scaling,
## but with NO regional averages and NO ratio printed, and panel letters re-shuffled with a fixed
## seed so the rater cannot match them to S11 Fig B. Key written to threshold-examples-blind-key.csv.
## Rebuilds the plug-in rate fields from the published common-minimum matrix, as analyze-systematic.R does.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript threshold-examples-blind.R
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
D <- readRDS("design-regions.rds"); ex <- fread("threshold-examples.csv")
mn <- readRDS(file.path(repo, "code-spatial-smoothing/controls/sweep/m-est-min-matrix.rds"))
keep <- !grepl("conflicted", mn$features); Mpub <- mn$counts[keep, ]; rownames(Mpub) <- mn$features[keep]
tc <- readRDS(file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
stopifnot(all(colnames(Mpub) == paste(tc$x, tc$y, sep = ".")))
R <- sweep(Mpub, 2, tc$total.count, "/")
xy <- data.table(i = (tc$x - 15025) / 50 + 1, j = (tc$y - 11525) / 50 + 1)
edges <- function(nm) { r <- D$regions[[nm]][, .(i, j)]; inr <- paste(r$i, r$j)
  rbind(r[!paste(i - 1, j) %in% inr, .(x = i - .5, xend = i - .5, y = j - .5, yend = j + .5)],
        r[!paste(i + 1, j) %in% inr, .(x = i + .5, xend = i + .5, y = j - .5, yend = j + .5)],
        r[!paste(i, j - 1) %in% inr, .(x = i - .5, xend = i + .5, y = j - .5, yend = j - .5)],
        r[!paste(i, j + 1) %in% inr, .(x = i - .5, xend = i + .5, y = j + .5, yend = j + .5)]) }
set.seed(20261029); ex <- ex[sample(.N)]; ex[, blind := LETTERS[seq_len(.N)]]
fwrite(ex[, .(blind, original = label, pair, r1, r2, ratio = round(ratio, 2))], "threshold-examples-blind-key.csv")
panels <- lapply(seq_len(nrow(ex)), function(k) { e <- ex[k]; v <- R[e$pair, ]; cap <- quantile(v, 0.99)
  d <- copy(xy)[, val := pmin(v, cap)]
  ggplot(d, aes(i, j)) + geom_raster(aes(fill = val)) +
    geom_segment(data = edges(e$r1), aes(x = x, xend = xend, y = y, yend = yend), color = "white", linewidth = 0.4) +
    geom_segment(data = edges(e$r2), aes(x = x, xend = xend, y = y, yend = yend), color = "red", linewidth = 0.4) +
    scale_fill_viridis_c(name = "rate") + coord_fixed(expand = FALSE) + theme_void(8) +
    labs(title = e$blind, subtitle = "white region vs red region: is the LEVEL of the field different?") })
ggsave("threshold-examples-blind.png", wrap_plots(panels, ncol = 2), width = 12, height = 13, dpi = 120)
cat("written\n")
