## Regenerate the Figure 3 and Figure 4 panel grids with a correct 1:1 aspect ratio.
##
## Why this script exists (reviewer comment R2.6): the published panels were drawn with
## `asp = qp.res.y / qp.res.x` (= 130/200 = 0.65). In R base graphics `asp` is the ratio of
## y-unit length to x-unit length, so 0.65 compresses y relative to x and stretches every
## bin horizontally by 1/0.65 = 1.54x. The bins are genuinely square (the binning steps
## x and y by the same `agg.size`, and the Stereo-seq spot pitch is isotropic), so the
## correct setting is `asp = 1`. That has now been corrected at source in the pipeline
## scripts; this script rebuilds only the two figures the reviewer flagged.
##
## Run:  distrobox enter emacs-r -- Rscript regenerate-figs-3-4.R
## Inputs:  data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-{poi,zip,zap}.rds
## Outputs: Fig3-raw-vs-normalized.png, Fig4-likelihood-comparison.png (300 dpi)

suppressMessages({library(fields); library(data.table)})

repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
pdir <- file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects")
odir <- file.path(repo, "Submission_PLOSOne/submission-2/figures-regenerated")
feat <- "Cd44"

get <- function(lik) readRDS(file.path(pdir, sprintf("%s-pred-%s.rds", feat, lik)))
poi <- get("poi"); zip <- get("zip"); zap <- get("zap")

d  <- as.data.table(poi$expect)
nx <- length(unique(d$x)); ny <- length(unique(d$y))
stopifnot(nx == 200, ny == 130)

## square bins: one x unit and one y unit are the same physical distance
ASP <- 1

panel <- function(dt, z, ttl, zlim = NULL) {
  fields::quilt.plot(dt$x, dt$y, z, nx = nx, ny = ny, asp = ASP,
                     main = ttl, xlab = "", ylab = "",
                     zlim = zlim, cex.main = 1.5)
}

## ---------------------------------------------------------------- Figure 3
## Raw vs count-normalized topology. Columns: feature counts, feature/total, total.
## Rows: observed data, model estimate, residual (data - model).
png(file.path(odir, "Fig3-raw-vs-normalized.png"),
    width = (nx/ny) * 6 * 3, height = 6 * 3 + 1, units = "in", res = 300)
par(mfrow = c(3, 3), mar = c(3, 3, 4, 6), oma = c(0, 0, 3, 0))

obs_raw  <- d$feat.count
obs_norm <- d$feat.count / d$total.count
tot      <- d$total.count
est_raw  <- d$mean
est_norm <- d$mean / d$total.count

panel(d, obs_raw,  sprintf("Observed %s counts", feat))
panel(d, obs_norm, sprintf("Observed %s / total counts", feat))
panel(d, tot,      "Total counts (nUMI)")

panel(d, est_raw,  sprintf("Modeled %s counts", feat))
panel(d, est_norm, sprintf("Modeled %s / total counts", feat))
panel(d, tot,      "Total counts (nUMI)")

panel(d, obs_raw  - est_raw,  "Residual, counts")
panel(d, obs_norm - est_norm, "Residual, normalized")
plot.new()  # deliberately blank: the total counts are an observed offset, not modeled,
            # so there is no residual for the third column

mtext(sprintf("Raw vs count-normalized topology (%s), 1:1 aspect ratio", feat),
      outer = TRUE, cex = 1.6, font = 2)
dev.off()

## ---------------------------------------------------------------- Figure 4
## Likelihood comparison, matching the published 4 x 4 layout.
## Rows: Poisson, ZIP, ZAP, raw data. Columns: P(non-zero), density (counts per nUMI),
## expected counts, residual. The bottom-right panel is the nUMI, as in the published
## caption, since the data row has no residual.
##
## Column 1 is recoverable for all three models; see export-panels.R for the derivation
## and its verification. In short: Poisson is exact from the expected count; ZIP's
## zero-inflation is a scalar hyperparameter recoverable as expect/(lambda*total.count);
## and ZAP's spatially varying presence field was saved all along.

p_nonzero <- function(o, which) {
  ex <- as.data.table(o$expect); lm_ <- as.data.table(o$lambda)
  switch(which,
    poi = 1 - exp(-ex$mean),
    zip = { mu <- lm_$mean * lm_$total.count
            median(ex$mean / mu) * (1 - exp(-mu)) },
    zap = as.data.table(o$presence)$mean)
}

png(file.path(odir, "Fig4-likelihood-comparison.png"),
    width = (nx/ny) * 6 * 4, height = 6 * 4 + 1, units = "in", res = 300)
par(mfrow = c(4, 4), mar = c(3, 3, 4, 6), oma = c(0, 0, 3, 0))

for (nm in c("Poisson", "ZIP", "ZAP")) {
  o   <- switch(nm, Poisson = poi, ZIP = zip, ZAP = zap)
  key <- switch(nm, Poisson = "poi", ZIP = "zip", ZAP = "zap")
  ex  <- as.data.table(o$expect); lm_ <- as.data.table(o$lambda)
  panel(ex, p_nonzero(o, key),       sprintf("%s: P(non-zero)", nm), zlim = c(0, 1))
  panel(ex, lm_$mean,                sprintf("%s: density (counts/nUMI)", nm))
  panel(ex, ex$mean,                 sprintf("%s: expected counts", nm))
  panel(ex, ex$feat.count - ex$mean, sprintf("%s: residual", nm))
}

## row 4: the raw data, for comparison against the three model rows above it
panel(d, as.numeric(d$feat.count > 0), "Data: observed presence", zlim = c(0, 1))
panel(d, d$feat.count / d$total.count, "Data: density (counts/nUMI)")
panel(d, d$feat.count,                 "Data: observed counts")
panel(d, d$total.count,                "Data: total counts (nUMI)")

mtext(sprintf("Likelihood model comparison (%s), 1:1 aspect ratio", feat),
      outer = TRUE, cex = 1.6, font = 2)
dev.off()

cat("wrote Fig3-raw-vs-normalized.png and Fig4-likelihood-comparison.png to\n  ", odir, "\n")
