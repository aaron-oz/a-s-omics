## Five plug-in examples for AOZ to rule on "looks different" (Phase 2). PLUG-IN ONLY.
## One instance drawn at random (seed 20261014) near each target ratio r = 1.15, 1.3, 1.5, 2, 3,
## from all (pair, frozen region comparison) instances with |r / target - 1| < 0.03.
## Maps: plug-in rate field (common minimum of posterior medians / total count), linear color
## scale capped at the field's 99th percentile, both regions outlined.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript threshold-examples.R
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
D <- readRDS("design-regions.rds")
P <- readRDS(file.path(Sys.getenv("TMPOUT", "."), "plugin-rate.rds"))   # written by plugin-ratios.R (about 300 MB; keep out of git)
out <- fread("plugin-ratios.csv")
set.seed(20261014)
targets <- c(1.15, 1.3, 1.5, 2, 3)
ex <- rbindlist(lapply(targets, function(t) { c <- out[abs(ratio / t - 1) < 0.03]; c[sample(.N, 1)][, target := t] }))
ex <- ex[sample(.N)]            # shuffle so the panels are not in ratio order
ex[, label := LETTERS[seq_len(.N)]]
print(ex[, .(label, pair, r1, r2, a = signif(a, 3), b = signif(b, 3), ratio = round(ratio, 2))])
fwrite(ex, "threshold-examples.csv")
xy <- data.table(x = as.numeric(sub("\\..*", "", P$key)), y = as.numeric(sub(".*\\.", "", P$key)))
xy[, `:=`(i = (x - 15025)/50 + 1, j = (y - 11525)/50 + 1)]
edges <- function(nm) {             # boundary segments of a region's bins
  r <- D$regions[[nm]][, .(i, j)]; inr <- paste(r$i, r$j)
  rbind(r[!paste(i - 1, j) %in% inr, .(x = i - .5, xend = i - .5, y = j - .5, yend = j + .5)],
        r[!paste(i + 1, j) %in% inr, .(x = i + .5, xend = i + .5, y = j - .5, yend = j + .5)],
        r[!paste(i, j - 1) %in% inr, .(x = i - .5, xend = i + .5, y = j - .5, yend = j - .5)],
        r[!paste(i, j + 1) %in% inr, .(x = i - .5, xend = i + .5, y = j + .5, yend = j + .5)])
}
panels <- lapply(seq_len(nrow(ex)), function(k) {
  e <- ex[k]; v <- P$R[e$pair, ]; cap <- quantile(v, 0.99)
  d <- copy(xy)[, val := pmin(v, cap)]
  ggplot(d, aes(i, j)) + geom_raster(aes(fill = val)) +
    geom_segment(data = edges(e$r1), aes(x = x, xend = xend, y = y, yend = yend), color = "white", linewidth = 0.4) +
    geom_segment(data = edges(e$r2), aes(x = x, xend = xend, y = y, yend = yend), color = "red", linewidth = 0.4) +
    scale_fill_viridis_c(name = "rate") + coord_fixed(expand = FALSE) + theme_void(8) +
    labs(title = sprintf("%s: %s", e$label, e$pair),
         subtitle = sprintf("white %s: avg %.3g   red %s: avg %.3g   ratio %.2f", e$r1, e$a, e$r2, e$b, e$ratio))
})
ggsave("threshold-examples.png", wrap_plots(panels, ncol = 2), width = 12, height = 13, dpi = 120)
