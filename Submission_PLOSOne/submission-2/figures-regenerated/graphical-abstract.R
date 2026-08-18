## Draft graphical abstract (R1.5).
##
## Reviewer 1 asked for something better than panels copy-pasted from the data section, so
## this is an argument rather than a montage: the headline result at left, the three controls
## that constrain what it means down the right, and the depth dependence along the bottom.
## Every number shown is measured. Nothing here is illustrative.
##
## Run: distrobox enter emacs-r -- Rscript graphical-abstract.R [width_in] [dpi]

suppressMessages(library(data.table))

a   <- commandArgs(trailingOnly = TRUE)
W   <- if (length(a) >= 1) as.numeric(a[1]) else 7.5
DPI <- if (length(a) >= 2) as.numeric(a[2]) else 400
H   <- 5.0

repo  <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
ctrl  <- file.path(repo, "data-outputs/mouse-embryo-raw/controls")
sweep <- file.path(repo, "code-spatial-smoothing/controls/sweep")
odir  <- file.path(repo, "Submission_PLOSOne/submission-2/figures-regenerated")

e <- new.env()
load(file.path(repo, "sam_sandbox/clustering vs features/m.est.min",
               "clustering vs features experiment",
               "meta.data.feature.depth.experiment.2025-07-12.Robj"), envir = e)
md <- get(ls(e)[1], envir = e)
cc <- grep("featureClusters$", colnames(md), value = TRUE)
depth_col <- function(n) cc[which(as.integer(sub("featureClusters$", "", cc)) == n)]

PAL <- unique(c('firebrick','steelblue','springgreen','purple','salmon','skyblue','navyblue',
  'violetred','#40376E','grey20','tomato','sandybrown','saddlebrown','royalblue','plum4',
  'lightgoldenrod','lawngreen','forestgreen','dimgray','deeppink','red2','paleturquoise1',
  'palevioletred','#E0CA3C','#989FCE','#261132','#F37748','#AFFC41','#347FC4','#0B5351',
  '#00A9A5','#4E8098','#D81159','#AB3428','#01BAEF','#F42C04','#610345','#044B7F','orchid4',
  'purple4','plum1','olivedrab2','slateblue','mediumvioletred','sienna','orange','seagreen',
  'lightseagreen','mediumpurple4','#F8F4A6','yellow','#A88FAC','#D2D68D','#CBD4C2'))

## Render as a raster rather than as points. The bins are a regular 200 x 130 grid, so
## plotting them as pch = 15 leaves visible gaps at any point size that does not also
## overlap, and the map reads as noise rather than as domains.
UX <- sort(unique(md$x)); UY <- sort(unique(md$y))
domain_map <- function(n, title = NULL, tcex = 0.62, mar = c(0.1, 0.1, 0.9, 0.1)) {
  lab <- as.integer(as.factor(md[[depth_col(n)]]))
  k <- max(lab)
  m <- matrix(NA_integer_, length(UX), length(UY))
  m[cbind(match(md$x, UX), match(md$y, UY))] <- lab
  par(mar = mar)
  image(UX, UY, m, col = PAL[(seq_len(k) - 1) %% length(PAL) + 1], zlim = c(1, k),
        asp = 1, axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  ## border at the data extent, not the plot region: asp = 1 centres the raster and a box()
  ## on the region leaves visible white bands above and below it
  dx <- diff(UX)[1] / 2; dy <- diff(UY)[1] / 2
  rect(min(UX) - dx, min(UY) - dy, max(UX) + dx, max(UY) + dy, lwd = 0.5, border = "grey35")
  if (!is.null(title)) mtext(title, side = 3, line = 0.1, cex = tcex, font = 2)
  invisible(k)
}

## the shape all three controls share: a null distribution with the observed value marked
verdict <- function(nulls, obs, headline, sub, xlab) {
  par(mar = c(1.55, 0.4, 1.95, 0.4))
  rng <- range(c(nulls, obs)); pad <- 0.12 * diff(rng)
  h <- hist(nulls, breaks = 20, plot = FALSE)
  plot(NA, xlim = c(rng[1] - pad, rng[2] + pad), ylim = c(0, max(h$counts) * 1.30),
       axes = FALSE, xlab = "", ylab = "")
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$counts,
       col = "grey84", border = "grey62", lwd = 0.4)
  abline(v = obs, col = "#C1272D", lwd = 2.3)
  axis(1, cex.axis = 0.54, tcl = -0.14, mgp = c(3, 0.14, 0), lwd = 0.6)
  mtext(xlab, side = 1, line = 0.68, cex = 0.46, col = "grey35")
  mtext(headline, side = 3, line = 0.80, cex = 0.64, font = 2)
  mtext(sub, side = 3, line = 0.18, cex = 0.48, col = "grey25")
  usr <- par("usr")
  adj <- if (obs < mean(usr[1:2])) c(-0.10, 0.5) else c(1.10, 0.5)
  text(obs, max(h$counts) * 1.19, "observed", col = "#C1272D", cex = 0.50, font = 2, adj = adj)
}

## ---- measured inputs
b <- rbindlist(lapply(list.files(file.path(ctrl, "control-B-klone"), "csv$", full.names = TRUE), fread))
cA <- rbindlist(lapply(list.files(file.path(sweep, "results-ctrlA"), "csv$", full.names = TRUE), fread))
cC <- readRDS(file.path(ctrl, "control-C-result.rds"))
sf <- unlist(lapply(c("results-seurat-d10","results-seurat-d50","results-seurat",
                      "results-seurat-350","results-seurat-d650","results-seurat-d1000"),
                    function(d) list.files(file.path(sweep, d), "csv$", full.names = TRUE)))
sw <- rbindlist(lapply(sf, fread))
dc <- sw[arm == "random", .(m = mean(elsa_d50), s = sd(elsa_d50)), by = depth][order(depth)]

## ------------------------------------------------------------------ compose
png(file.path(odir, "graphical-abstract.png"), width = W, height = H, units = "in", res = DPI)
L <- rbind(rep(1, 12),
           c(rep(2, 7), rep(3, 5)),
           c(rep(2, 7), rep(4, 5)),
           c(rep(2, 7), rep(5, 5)),
           c(rep(6, 3), rep(7, 3), rep(8, 3), rep(9, 3)))
layout(L, heights = c(0.30, 1, 1, 1, 1.15))
par(omi = c(0.05, 0.05, 0.05, 0.05))

## 1 title
par(mar = c(0, 0, 0, 0)); plot.new()
text(0.5, 0.68, "Signalling-field domains recover tissue architecture", cex = 1.12, font = 2)
text(0.5, 0.22, paste("Continuous ligand-receptor fields inferred from raw counts, clustered without",
                      "dimensionality reduction"), cex = 0.62, col = "grey25")

## 2 headline domain map
nd <- domain_map(1481, NULL, mar = c(0.1, 0.1, 1.0, 0.1))
mtext(sprintf("1,477 interaction fields  →  %d spatial domains", nd),
      side = 3, line = 0.20, cex = 0.76, font = 2)

## 3-5 the three controls
verdict(b[perm_id > 0, elsa_d50], b[perm_id == 0, elsa_d50],
        "Not made by the smoother",
        sprintf("coordinates permuted, all fields refit: %d domains collapse to %d-%d",
                b[perm_id == 0, n_clusters], min(b[perm_id > 0, n_clusters]),
                max(b[perm_id > 0, n_clusters])),
        "mean ELSA (d = 50), lower = more coherent")
verdict(cC$null, cC$ari, "Matches independent anatomy",
        sprintf("ARI %.3f against annotation the pipeline never saw", cC$ari),
        "adjusted Rand index vs 499 label permutations")
verdict(cA[tag != "observed", elsa_d50], cA[tag == "observed", elsa_d50],
        "Robust to how genes are paired",
        "99 scrambled pairings, no real pair retained",
        "mean ELSA (d = 50)")

## 6-8 depth strip, 9 the deconfounded curve
for (n in c(10, 150, 1481)) {
  k <- domain_map(n, NULL, mar = c(0.35, 0.1, 0.85, 0.1))
  mtext(sprintf("%s fields", format(if (n == 1481) 1477 else n, big.mark = ",")),
        side = 3, line = 0.05, cex = 0.58, font = 2)
}
par(mar = c(1.9, 2.0, 0.85, 0.4))
plot(dc$depth, dc$m, type = "n", ylim = c(0, max(dc$m + dc$s)), axes = FALSE,
     xlab = "", ylab = "", log = "x")
arrows(dc$depth, dc$m - dc$s, dc$depth, dc$m + dc$s, code = 3, angle = 90,
       length = 0.015, col = "grey65", lwd = 0.8)
lines(dc$depth, dc$m, col = "#22577A", lwd = 1.8)
points(dc$depth, dc$m, pch = 19, cex = 0.6, col = "#22577A")
axis(1, at = c(10, 100, 1000), labels = c("10", "100", "1,000"),
     cex.axis = 0.55, tcl = -0.14, mgp = c(3, 0.12, 0), lwd = 0.6)
axis(2, cex.axis = 0.55, tcl = -0.14, mgp = c(3, 0.30, 0), lwd = 0.6, las = 1)
mtext("number of fields (random subsets)", side = 1, line = 0.85, cex = 0.48, col = "grey35")
mtext("mean ELSA", side = 2, line = 1.15, cex = 0.48, col = "grey35")
mtext("Needs hundreds, not a handful", side = 3, line = 0.05, cex = 0.58, font = 2)
dev.off()
cat(sprintf("wrote graphical-abstract.png  %.1f x %.1f in @ %d dpi, %d domains\n", W, H, DPI, nd))
