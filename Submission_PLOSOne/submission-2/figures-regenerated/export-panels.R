## Export individual sub-panels for Figures 3 and 4 at a correct 1:1 aspect ratio,
## as separate files for placement in Adobe Illustrator.
##
## Each panel is the map plus its colour bar, with NO title and minimal margin, so it can be
## dropped straight into the existing figure layout. Panel identity is encoded in the
## filename and tabulated in README-panels.md.
##
## Aspect: the field of view is 200 x 130 bins and the bins are square, so asp = 1 is correct.
## The published panels used asp = qp.res.y/qp.res.x = 0.65, which stretched every bin
## horizontally by 1.54. See ../../../ for the source fix.
##
## Run: distrobox enter emacs-r -- Rscript export-panels.R

suppressMessages({library(fields); library(data.table)})

repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
pdir <- file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects")
odir <- file.path(repo, "Submission_PLOSOne/submission-2/figures-regenerated/panels")
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
feat <- "Cd44"

get <- function(lik) readRDS(file.path(pdir, sprintf("%s-pred-%s.rds", feat, lik)))
poi <- get("poi"); zip <- get("zip"); zap <- get("zap")

d  <- as.data.table(poi$expect)
nx <- length(unique(d$x)); ny <- length(unique(d$y))
stopifnot(nx == 200, ny == 130)

## device sized so the map area is close to the 200:130 domain aspect once the
## colour bar and margins are accounted for
emit <- function(fname, z, dt = d) {
  png(file.path(odir, paste0(fname, ".png")),
      width = 8, height = 5.0, units = "in", res = 300)
  par(mar = c(2.2, 2.2, 0.6, 1.0))
  fields::quilt.plot(dt$x, dt$y, z, nx = nx, ny = ny, asp = 1,
                     xlab = "", ylab = "", cex.axis = 0.8)
  dev.off()
  invisible(fname)
}

## ---- Figure 3: raw vs count-normalized -------------------------------------
## cols: feature counts | feature/total | total ;  rows: observed | modeled | residual
emit("Fig3_r1c1_observed-counts",        d$feat.count)
emit("Fig3_r1c2_observed-normalized",    d$feat.count / d$total.count)
emit("Fig3_r1c3_total-counts-nUMI",      d$total.count)
emit("Fig3_r2c1_modeled-counts",         d$mean)
emit("Fig3_r2c2_modeled-normalized",     d$mean / d$total.count)
emit("Fig3_r2c3_total-counts-nUMI",      d$total.count)
emit("Fig3_r3c1_residual-counts",        d$feat.count - d$mean)
emit("Fig3_r3c2_residual-normalized",    d$feat.count/d$total.count - d$mean/d$total.count)
## r3c3 is deliberately empty: the total counts are an observed offset, not modeled.

## ---- Figure 4: likelihood comparison ---------------------------------------
## rows: Poisson | ZIP | ZAP | raw ;  cols: P(non-zero) | density | counts | residual
for (nm in c("poi","zip","zap")) {
  o  <- switch(nm, poi = poi, zip = zip, zap = zap)
  ex <- as.data.table(o$expect); lm_ <- as.data.table(o$lambda)
  tag <- switch(nm, poi = "r1_Poisson", zip = "r2_ZIP", zap = "r3_ZAP")
  emit(sprintf("Fig4_%s_c2_density", tag),  lm_$mean)
  emit(sprintf("Fig4_%s_c3_counts",  tag),  ex$mean)
  emit(sprintf("Fig4_%s_c4_residual",tag),  ex$feat.count - ex$mean)
}
## Column 1, probability of a non-zero value. Exact for the Poisson model, where
## P(Y>0) = 1 - exp(-mu) and mu is the expected count. NOT reconstructable for ZIP and
## ZAP: the saved prediction objects store lambda, expect and obs_prob only, so the
## spatially varying zero probability p_i(s) is not available. Those two panels must come
## from a rerun that saves p_i(s), or be dropped from the figure.
emit("Fig4_r1_Poisson_c1_prob-nonzero", 1 - exp(-d$mean))

## row 4: the raw data
emit("Fig4_r4_raw_c1_presence",     as.numeric(d$feat.count > 0))
emit("Fig4_r4_raw_c2_density",      d$feat.count / d$total.count)
emit("Fig4_r4_raw_c3_counts",       d$feat.count)
emit("Fig4_r4_raw_c4_total-nUMI",   d$total.count)

cat("wrote", length(list.files(odir, pattern="[.]png$")), "panels to", odir, "\n")
