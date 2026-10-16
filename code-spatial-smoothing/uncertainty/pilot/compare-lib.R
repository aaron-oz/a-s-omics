## Comparison functions for the uncertainty pivot (Phase 3, step 4). Definitions follow
## DESIGN.org "Quantities compared" and "Inference". All on the rate scale; contrasts are log ratios.
## Draw matrices are bins x draws (26,000 x D), rows in the saved bin order (x, then y).
suppressMessages({library(data.table); library(excursions)})
REPO <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
DDIR <- "/var/home/aoz/data/a-s-omics-uncertainty/draws"
UDIR <- file.path(REPO, "code-spatial-smoothing/uncertainty")

load_draws <- function(g) readRDS(file.path(DDIR, paste0(g, "-lambda-draws.rds")))$draws
gene_case <- function(G) paste0(substr(G, 1, 1), tolower(substring(G, 2)))   # CALR -> Calr

## Posterior draws of a signaling field: M^(d)(s) = min(lambda_L^(d)(s), lambda_R^(d)(s))
sig_draws <- function(pair, cache) {
  lr <- strsplit(pair, "-")[[1]]
  pmin(cache[[gene_case(lr[1])]], cache[[gene_case(lr[2])]])
}
## Shape normalization within each draw: S^(d)(s) = M^(d)(s) / mean over all bins of M^(d)
shape_norm <- function(M) sweep(M, 2, colMeans(M), "/")

## Areal-average contrast between bin sets k1, k2: log A_k1 - log A_k2, per draw
area_contrast <- function(M, k1, k2) log(colMeans(M[k1, , drop = FALSE])) - log(colMeans(M[k2, , drop = FALSE]))

summ_contrast <- function(c) {
  q <- quantile(c, c(.025, .5, .975))
  list(post_ratio = exp(q[[2]]), lo = exp(q[[1]]), hi = exp(q[[3]]),
       p_gt1 = mean(c > 0), credible = (q[[1]] > 0) || (q[[3]] < 0))
}

## Along a slice: same-field deviation from the slice average, per draw
slice_dev <- function(M, k) { X <- log(M[k, , drop = FALSE]); sweep(X, 2, log(colMeans(M[k, , drop = FALSE])), "-") }

## Joint and pointwise inference on a matrix of contrast draws (locations x draws)
joint_summary <- function(X, alpha = 0.05) {
  pos <- excursions.mc(X, alpha = alpha, u = 0, type = ">")
  neg <- excursions.mc(X, alpha = alpha, u = 0, type = "<")
  pw  <- t(apply(X, 1, quantile, c(alpha / 2, .5, 1 - alpha / 2)))
  list(pos = as.logical(pos$E), neg = as.logical(neg$E), Fpos = pos$F, Fneg = neg$F,
       pw_lo = pw[, 1], med = pw[, 2], pw_hi = pw[, 3],
       ## pointwise calls at the same one-sided level as the excursion sets, so E is a subset of them
       pw_pos = rowMeans(X > 0) >= 1 - alpha, pw_neg = rowMeans(X < 0) >= 1 - alpha)
}
