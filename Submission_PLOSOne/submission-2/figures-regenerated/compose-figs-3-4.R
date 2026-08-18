## Build Figures 3 and 4 as complete, publication-sized composites from a script.
##
## Replaces hand assembly in Illustrator. The published composites had no layout script, which
## the manuscript had to admit under R2.9; this file is that script. It makes the R1.3 sizing a
## parameter rather than a manual pass, and it inherits the R2.6 fix (asp = 1) by construction.
##
## Two things buy the panels their size, which is what R1.3 actually asked for. Every column
## uses one common colour scale, so a single shared colour bar per column replaces sixteen
## per-panel bars and returns that width to the maps. And the figure height is solved from the
## panel geometry so the plot region matches the 200:130 data aspect exactly, leaving no dead
## space above or below each map.
##
## Run: distrobox enter emacs-r -- Rscript compose-figs-3-4.R [width_in] [dpi]

suppressMessages({library(fields); library(data.table)})

a   <- commandArgs(trailingOnly = TRUE)
W   <- if (length(a) >= 1) as.numeric(a[1]) else 7.5   # PLOS full-width maximum
DPI <- if (length(a) >= 2) as.numeric(a[2]) else 400

repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
pdir <- file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects")
odir <- file.path(repo, "Submission_PLOSOne/submission-2/figures-regenerated")
feat <- "Cd44"

get <- function(lik) readRDS(file.path(pdir, sprintf("%s-pred-%s.rds", feat, lik)))
poi <- get("poi"); zip <- get("zip"); zap <- get("zap")
d  <- as.data.table(poi$expect)
nx <- length(unique(d$x)); ny <- length(unique(d$y))
stopifnot(nx == 200, ny == 130)

p_nonzero <- function(o, which) {           # derivations in export-panels.R
  ex <- as.data.table(o$expect); lm_ <- as.data.table(o$lambda)
  switch(which,
    poi = 1 - exp(-ex$mean),
    zip = { mu <- lm_$mean * lm_$total.count; median(ex$mean / mu) * (1 - exp(-mu)) },
    zap = as.data.table(o$presence)$mean)
}

CBAR_IN <- 0.40                    # height of the shared colour-bar strip, inches
OMI <- c(0.06 + CBAR_IN, 0.46, 0.34, 0.06)   # bottom (incl. bar strip), left, top, right
MAI <- c(0.03, 0.03, 0.03, 0.03)
COLS <- tim.colors(256)

## Solve the figure height so each plot region is exactly the data aspect: no dead space.
solve_h <- function(ncol, nrow) {
  plot_w <- (W - OMI[2] - OMI[4]) / ncol - MAI[2] - MAI[4]
  plot_h <- plot_w * ny / nx
  list(H = nrow * (plot_h + MAI[1] + MAI[3]) + OMI[1] + OMI[3], pw = plot_w)
}

LAB <- 1
panel <- function(z, zlim = NULL, blank = FALSE) {
  if (blank) { plot.new(); return(invisible()) }
  fields::quilt.plot(d$x, d$y, z, nx = nx, ny = ny, asp = 1, zlim = zlim,
                     xlab = "", ylab = "", axes = FALSE, add.legend = FALSE, col = COLS)
  box(lwd = 0.5, col = "grey35")
  usr <- par("usr")
  text(usr[1] + 0.028 * diff(usr[1:2]), usr[4] - 0.055 * diff(usr[3:4]),
       LETTERS[LAB], adj = c(0, 1), font = 2, cex = 1.0)
  LAB <<- LAB + 1
}

## `smallplot` in fields::image.plot is in DEVICE coordinates, not cell coordinates, so the
## bars must be positioned explicitly rather than by advancing through a layout row.
cbar <- function(zlim, i, ncol, H) {
  cw <- (W - OMI[2] - OMI[4]) / ncol
  x0 <- (OMI[2] + (i - 1) * cw + 0.16 * cw) / W
  x1 <- (OMI[2] + (i - 1) * cw + 0.84 * cw) / W
  y0 <- (0.06 + 0.56 * CBAR_IN) / H
  y1 <- (0.06 + 0.80 * CBAR_IN) / H
  ## smallplot is in DEVICE coordinates, so reset to a full-device figure region first.
  ## image.plot(legend.only) also leaves par pointing at the legend it just drew, so this
  ## reset is needed before every bar, not only the first.
  par(new = TRUE, fig = c(0, 1, 0, 1), omi = c(0, 0, 0, 0), mai = c(0, 0, 0, 0))
  plot.new()
  image.plot(zlim = zlim, legend.only = TRUE, horizontal = TRUE, col = COLS,
             smallplot = c(x0, x1, y0, y1),
             axis.args = list(cex.axis = 0.60, mgp = c(3, 0.20, 0), tcl = -0.14))
}

## Units live in the column headers rather than under the colour bars. Labelling the bars
## separately means placing text against image.plot's own coordinate system, which is not
## reliably recoverable, and the header is the less cluttered place for them regardless.
colhead <- function(labs) for (i in seq_along(labs))
  mtext(labs[i], side = 3, outer = TRUE, cex = 0.80, font = 2, line = 0.35,
        at = (OMI[2] + (i - 0.5) * (W - OMI[2] - OMI[4]) / length(labs)) / W)
rowlab <- function(labs, nrow, H) {
  usable <- H - OMI[1] - OMI[3]
  for (i in seq_along(labs))
    mtext(labs[i], side = 2, outer = TRUE, cex = 0.80, font = 2, line = 0.5,
          at = (OMI[1] + usable - (i - 0.5) * usable / nrow) / H)
}

## ------------------------------------------------------------------ Figure 3
g <- solve_h(3, 3); H3 <- g$H
png(file.path(odir, "Fig3-composite.png"), width = W, height = H3, units = "in", res = DPI)
par(omi = OMI, mai = MAI, mfrow = c(3, 3)); LAB <- 1

obs_raw <- d$feat.count; obs_norm <- d$feat.count / d$total.count; tot <- d$total.count
est_raw <- d$mean;       est_norm <- d$mean / d$total.count
zl_raw  <- range(c(obs_raw, est_raw)); zl_norm <- range(c(obs_norm, est_norm))
zl_tot  <- range(tot)
zl_rraw <- range(obs_raw - est_raw);   zl_rnorm <- range(obs_norm - est_norm)

panel(obs_raw, zl_raw); panel(obs_norm, zl_norm); panel(tot, zl_tot)
panel(est_raw, zl_raw); panel(est_norm, zl_norm); panel(tot, zl_tot)
panel(obs_raw - est_raw, zl_rraw); panel(obs_norm - est_norm, zl_rnorm); panel(blank = TRUE)
colhead(c(sprintf("%s counts", feat), sprintf("%s / total counts", feat), "Total counts (nUMI)"))
rowlab(c("Observed", "Modeled", "Residual"), 3, H3)
cbar(zl_raw, 1, 3, H3); cbar(zl_norm, 2, 3, H3); cbar(zl_tot, 3, 3, H3)
dev.off()

## ------------------------------------------------------------------ Figure 4
g4 <- solve_h(4, 4); H4 <- g4$H
png(file.path(odir, "Fig4-composite.png"), width = W, height = H4, units = "in", res = DPI)
par(omi = OMI, mai = MAI, mfrow = c(4, 4)); LAB <- 1

dens <- lapply(list(poi, zip, zap), function(o) as.data.table(o$lambda)$mean)
cnts <- lapply(list(poi, zip, zap), function(o) as.data.table(o$expect)$mean)
zl_dens <- range(c(unlist(dens), d$feat.count / d$total.count))
zl_cnt  <- range(c(unlist(cnts), d$feat.count))
zl_res  <- range(unlist(lapply(cnts, function(m) d$feat.count - m)))

for (i in 1:3) {
  key <- c("poi", "zip", "zap")[i]; o <- list(poi, zip, zap)[[i]]
  panel(p_nonzero(o, key), c(0, 1)); panel(dens[[i]], zl_dens)
  panel(cnts[[i]], zl_cnt);          panel(d$feat.count - cnts[[i]], zl_res)
}
panel(as.numeric(d$feat.count > 0), c(0, 1)); panel(d$feat.count / d$total.count, zl_dens)
panel(d$feat.count, zl_cnt)
## the bottom-right panel is the nUMI, as in the published figure. It is on a different scale
## from the counts column above it, so it carries its own bar rather than sharing one.
fields::quilt.plot(d$x, d$y, d$total.count, nx = nx, ny = ny, asp = 1, xlab = "", ylab = "",
                   axes = FALSE, col = COLS, legend.width = 0.7, legend.mar = 3.0,
                   legend.cex = 0.5, axis.args = list(cex.axis = 0.5))
box(lwd = 0.5, col = "grey35")
usr <- par("usr"); text(usr[1] + 0.028 * diff(usr[1:2]), usr[4] - 0.055 * diff(usr[3:4]),
                        LETTERS[LAB], adj = c(0, 1), font = 2, cex = 1.0)
colhead(c("P(non-zero)", "Density (counts/nUMI)", "Counts", "Residual (obs - model)"))
rowlab(c("Poisson", "ZIP", "ZAP", "Observed"), 4, H4)
cbar(c(0, 1), 1, 4, H4);   cbar(zl_dens, 2, 4, H4)
cbar(zl_cnt, 3, 4, H4);    cbar(zl_res, 4, 4, H4)
dev.off()

cat(sprintf("Fig3-composite.png  %.2f x %.2f in @ %d dpi, panel width %.2f in\n", W, H3, DPI, g$pw))
cat(sprintf("Fig4-composite.png  %.2f x %.2f in @ %d dpi, panel width %.2f in\n", W, H4, DPI, g4$pw))
