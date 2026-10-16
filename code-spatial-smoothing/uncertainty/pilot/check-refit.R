## Reproduction check (Phase 3, step 3): does each refit reproduce the SAVED per-bin posterior
## medians of lambda (data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/<Gene>-pred-poi.rds)?
## Per bin: dlog = log(refit median) - log(saved median). Under pure Monte Carlo error, with
## both medians from 1000 draws and a roughly lognormal marginal with log-scale sd s (from the
## saved 2.5% and 97.5% quantiles), sd(dlog) is about sqrt(2) * 1.2533 * s / sqrt(1000).
## z = dlog / that sd should look like N(0, 1). Pass rule (fixed before running): |mean z| < 0.2,
## sd(z) < 1.3, and |median dlog| < 0.02. Also reports the refit vs saved width of the 95% interval.
## Usage: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript check-refit.R Cd44 [...]
suppressMessages(library(data.table))
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
pdir <- file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects")
ddir <- "/var/home/aoz/data/a-s-omics-uncertainty/draws"
res <- rbindlist(lapply(commandArgs(trailingOnly = TRUE), function(g) {
  d <- readRDS(file.path(ddir, paste0(g, "-lambda-draws.rds")))
  sv <- readRDS(file.path(pdir, paste0(g, "-pred-poi.rds")))$lambda
  stopifnot(all(sv$x == d$xy[, 1]), all(sv$y == d$xy[, 2]))
  q <- t(apply(d$draws, 1, quantile, c(.025, .5, .975)))
  s <- (log(sv$q0.975) - log(sv$q0.025)) / (2 * qnorm(.975))
  dlog <- log(q[, 2]) - log(sv$median)
  z <- dlog / (sqrt(2) * 1.2533 * s / sqrt(ncol(d$draws)))
  wr <- log(q[, 3] / q[, 1]) / log(sv$q0.975 / sv$q0.025)
  data.table(gene = g, med_dlog = median(dlog), q99_abs_dlog = quantile(abs(dlog), .99),
             mean_z = mean(z), sd_z = sd(z), width_ratio_med = median(wr),
             pass = abs(mean(z)) < 0.2 && sd(z) < 1.3 && abs(median(dlog)) < 0.02,
             range_med = d$hyper[1, "mean"], sigma_med = d$hyper[2, "mean"], secs = round(d$seconds))
}))
print(res, digits = 3)
fwrite(res, "check-refit.csv", append = file.exists("check-refit.csv"))
