## REVIEW-2026-10-29.org item P1: does the lognormal approximation behind the plug-in-bias section
## (jensen-bias.R: each marginal posterior of expected count taken as lognormal with mu = log q0.5,
## sigma = (log q0.975 - log q0.025) / (2 * 1.96); P(M > t) = P(L > t) P(R > t)) reproduce the
## posterior median of the interaction field computed from the joint draws? Checked on the seven
## showcase pairs, rate scale (expected count / total count), all 26,000 bins.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript lognormal-check.R
source("compare-lib.R")
pdir <- file.path(REPO, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects")
showcase <- c("CALR-ITGA2B", "NID1-ITGB1", "ADAM12-ITGB1", "WNT4-FZD6", "RSPO1-LGR6", "DHH-PTCH1", "VEGFA-KDR")
tc <- readRDS(file.path(pdir, "Cd44-pred-poi.rds"))$lambda; N <- tc$total.count
marg <- function(g) { ex <- readRDS(file.path(pdir, paste0(g, "-pred-poi.rds")))$expect
  stopifnot(all(ex$x == tc$x), all(ex$y == tc$y))
  list(mu = log(ex$median), sg = (log(ex$q0.975) - log(ex$q0.025)) / (2 * qnorm(.975))) }
## median of M under the approximation: solve P(L > t) P(R > t) = 1/2 by bisection on log t
approx_med <- function(a, b) { lo <- pmin(a$mu, b$mu) - 12; hi <- pmin(a$mu, b$mu) + 1
  for (k in 1:60) { m <- (lo + hi) / 2; p <- pnorm(m, a$mu, a$sg, lower.tail = FALSE) * pnorm(m, b$mu, b$sg, lower.tail = FALSE)
    up <- p > .5; lo[up] <- m[up]; hi[!up] <- m[!up] }
  exp((lo + hi) / 2) }
out <- rbindlist(lapply(showcase, function(f) { lr <- strsplit(f, "-")[[1]]; gl <- gene_case(lr[1]); gr <- gene_case(lr[2])
  cache <- setNames(list(load_draws(gl), load_draws(gr)), c(gl, gr))
  dm <- apply(sig_draws(f, cache), 1, median)                       # rate scale, from joint draws
  am <- approx_med(marg(gl), marg(gr)) / N                            # expected-count scale -> rate
  pl <- pmin(exp(marg(gl)$mu), exp(marg(gr)$mu)) / N                  # published plug-in field, rate scale
  d <- abs(log(dm) - log(am))
  data.table(pair = f, med_absdiff = median(d), p90_absdiff = quantile(d, .9), max_absdiff = max(d),
             pearson = cor(dm, am), median_over_draws_pct = 100 * median(pl / dm - 1), median_over_approx_pct = 100 * median(pl / am - 1),
             share_gt10_draws = mean(pl / dm - 1 > .1), share_gt10_approx = mean(pl / am - 1 > .1)) }))
fwrite(out, "results-lognormal-check.csv"); print(out, digits = 3)
