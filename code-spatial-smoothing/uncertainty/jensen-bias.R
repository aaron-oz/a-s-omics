## Plug-in bias of the common-minimum interaction fields.
##
## The published interaction field for a pair is m_plug(s) = min(med L(s), med R(s)), the minimum
## of the two posterior MEDIANS of expected count (03-array-post-process.R:167-168, 189). The
## coherent target is the posterior of the interaction itself, M(s) = min(L(s), R(s)).
## L and R come from separate model fits, so their posteriors are independent and, at each bin,
## M's distribution follows from the two marginals alone:
##     P(M > t) = P(L > t) P(R > t).
## Two summaries of M are computed:
##   med_M : posterior median of M. Exact inequality: m_plug >= med_M always, because with
##           m* = min(med L, med R), P(M > m*) = P(L > m*) P(R > m*) <= 1/2.
##   mean_M: posterior mean of M. min(E L, E R) >= E M by Jensen (min is concave); the published
##           field uses medians, not means, so its direction against mean_M is not guaranteed.
##
## Posterior draws were not saved, only summaries (mean, sd, q0.025, q0.5, q0.975 per bin).
## Each marginal is taken as lognormal, matching the model's Gaussian latent field on the log
## scale: mu = log(q0.5), sigma = (log q0.975 - log q0.025) / (2 * 1.959964). This assumption
## is checked against the saved posterior mean and sd before anything else is reported.
##
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript jensen-bias.R
## Output: jensen-bias-pairs.csv (one row per pair), jensen-bias-bins.rds (per-bin summaries),
##         jensen-bias.log via the caller.

suppressMessages(library(data.table))
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
pdir <- file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects")
sw   <- file.path(repo, "code-spatial-smoothing/controls/sweep")
odir <- file.path(repo, "code-spatial-smoothing/uncertainty")

## ---- pairs: the 1477 mechanisms of the published clustering, duplicates dropped
mn  <- readRDS(file.path(sw, "m-est-min-matrix.rds"))
keepm <- !grepl("conflicted", mn$features)
mech <- mn$features[keepm]; Mpub <- mn$counts[keepm, ]
e <- new.env(); load(file.path(sw, "fantom-hierarchy.Robj"), envir = e); fan <- get("output", e)
fan$key <- gsub("_", "-", fan$MECHANISM)
fi <- match(mech, fan$key)
stopifnot(!anyNA(fi))
Lg <- fan$LIGAND[fi]; Rg <- fan$RECEPTOR[fi]
files <- list.files(pdir, "-pred-poi\\.rds$")
gname <- sub("-pred-poi\\.rds$", "", files)
gmap <- setNames(gname, toupper(gname))
need <- unique(c(Lg, Rg)); stopifnot(all(need %in% names(gmap)))
cat(sprintf("pairs: %d   genes needed: %d\n", length(mech), length(need)))

## ---- read the needed marginals
NB <- ncol(Mpub)
MU <- SG <- MEAN <- SD <- matrix(NA_real_, NB, length(need), dimnames = list(NULL, need))
xy <- NULL
t0 <- Sys.time()
for (k in seq_along(need)) {
  ex <- readRDS(file.path(pdir, paste0(gmap[[need[k]]], "-pred-poi.rds")))$expect
  if (is.null(xy)) xy <- ex[, .(x, y)] else stopifnot(all(ex$x == xy$x), all(ex$y == xy$y))
  MU[, k]   <- log(ex$median)
  SG[, k]   <- (log(ex$q0.975) - log(ex$q0.025)) / (2 * qnorm(0.975))
  MEAN[, k] <- ex$mean; SD[, k] <- ex$sd
}
cat(sprintf("read %d genes in %.0f s\n", length(need), as.numeric(Sys.time() - t0, units = "secs")))

## ---- check 1: the published field is min of medians, with the same bin order
li <- match(Lg, need); ri <- match(Rg, need)
plug_check <- vapply(seq_along(mech)[1:50], function(i)
  max(abs(pmin(exp(MU[, li[i]]), exp(MU[, ri[i]])) - Mpub[i, ])), 0)
cat(sprintf("check 1, max |min(medians) - published| over 50 pairs: %.2e\n", max(plug_check)))

## ---- check 2: is lognormal adequate? compare implied to saved mean and sd
mean_ln <- exp(MU + SG^2 / 2)
sd_ln   <- mean_ln * sqrt(exp(SG^2) - 1)
rm_ <- as.vector(mean_ln / MEAN - 1); rs_ <- as.vector(sd_ln / SD - 1)
q <- function(v) quantile(v, c(.01, .05, .5, .95, .99), na.rm = TRUE)
cat("check 2, lognormal-implied / saved posterior MEAN - 1, quantiles over all gene-bins:\n"); print(round(q(rm_), 4))
cat("check 2, lognormal-implied / saved posterior SD - 1:\n"); print(round(q(rs_), 4))
rm(mean_ln, sd_ln, rm_, rs_); invisible(gc())

## ---- per-pair computation
GH <- local({            # 30-node Gauss-Hermite for E over a normal
  n <- 30; i <- seq_len(n - 1); b <- sqrt(i / 2)
  J <- matrix(0, n, n); J[cbind(i, i + 1)] <- b; J[cbind(i + 1, i)] <- b
  ev <- eigen(J, symmetric = TRUE)
  list(x = ev$values, w = ev$vectors[1, ]^2 * sqrt(pi))
})
## E[min(l, R)] for fixed l and lognormal R(mu, s)
emin_fixed <- function(l, mu, s) {
  z <- (log(l) - mu) / s
  exp(mu + s^2 / 2) * pnorm(z - s) + l * pnorm(z, lower.tail = FALSE)
}
pair_stats <- function(i) {
  mL <- MU[, li[i]]; sL <- SG[, li[i]]; mR <- MU[, ri[i]]; sR <- SG[, ri[i]]
  plug <- exp(pmin(mL, mR))
  ## posterior median of M: solve P(L>t)P(R>t) = 1/2 on log t, bracketed below plug
  lo <- pmin(mL, mR) - 8 * pmax(sL, sR); hi <- pmin(mL, mR)
  for (it in 1:45) {
    mid <- (lo + hi) / 2
    f <- pnorm((mid - mL) / sL, lower.tail = FALSE) * pnorm((mid - mR) / sR, lower.tail = FALSE)
    up <- f > 0.5; lo[up] <- mid[up]; hi[!up] <- mid[!up]
  }
  medM <- exp((lo + hi) / 2)
  ## posterior mean of M by Gauss-Hermite over L
  meanM <- 0
  for (k in seq_along(GH$x)) {
    l <- exp(mL + sqrt(2) * sL * GH$x[k])
    meanM <- meanM + GH$w[k] * emin_fixed(l, mR, sR)
  }
  meanM <- meanM / sqrt(pi)
  plug_mean <- pmin(exp(mL + sL^2 / 2), exp(mR + sR^2 / 2))   # min of the two posterior means
  zgap <- abs(mL - mR) / sqrt(sL^2 + sR^2)                     # standardized gap between partners
  list(plug = plug, medM = medM, meanM = meanM, plug_mean = plug_mean, zgap = zgap,
       slim = ifelse(mL <= mR, sL, sR))
}

bias_med_bins <- matrix(NA_real_, NB, length(mech))   # relative bias of published vs median of M
rows <- vector("list", length(mech))
set.seed(1); samp_pairs <- sample(seq_along(mech), 200)
scatter <- list()
t0 <- Sys.time()
for (i in seq_along(mech)) {
  r <- pair_stats(i)
  rb_med  <- r$plug / r$medM - 1
  rb_mean <- r$plug / r$meanM - 1
  jen     <- r$plug_mean / r$meanM - 1
  bias_med_bins[, i] <- rb_med
  rows[[i]] <- data.table(
    mech = mech[i],
    relbias_med_median = median(rb_med), relbias_med_p95 = quantile(rb_med, .95),
    frac_med_gt10 = mean(rb_med > 0.10),
    relbias_mean_median = median(rb_mean), frac_mean_neg = mean(rb_mean < 0),
    jensen_means_median = median(jen),
    cor_pattern_med  = cor(r$plug, r$medM),  cor_pattern_mean = cor(r$plug, r$meanM),
    spear_pattern_med = cor(r$plug, r$medM, method = "spearman"),
    cor_bias_zgap = cor(rb_med, r$zgap, method = "spearman"),
    cor_bias_slim = cor(rb_med, r$slim, method = "spearman"))
  if (i %in% samp_pairs) scatter[[length(scatter) + 1]] <-
    data.table(mech = mech[i], rb_med = rb_med, zgap = r$zgap, slim = r$slim)[sample(.N, 500)]
  if (i %% 250 == 0) cat(sprintf("  %d pairs, %.0f s\n", i, as.numeric(Sys.time() - t0, units = "secs")))
}
P <- rbindlist(rows)
fwrite(P, file.path(odir, "jensen-bias-pairs.csv"))
bins <- data.table(xy, median_relbias = apply(bias_med_bins, 1, median),
                   p90_relbias = apply(bias_med_bins, 1, quantile, .9))
saveRDS(list(bins = bins, scatter = rbindlist(scatter)), file.path(odir, "jensen-bias-bins.rds"))

cat("\n=== RESULTS ===\n")
cat(sprintf("min over all pairs and bins of published/median(M) - 1: %.3e (theory: >= 0)\n", min(bias_med_bins)))
cat("published vs posterior MEDIAN of M, relative bias, per-pair median over bins:\n"); print(round(q(P$relbias_med_median), 4))
cat("  pooled over all pair-bins:\n"); print(round(q(as.vector(bias_med_bins)), 4))
cat(sprintf("  fraction of pair-bins with relative bias > 10%%: %.3f; > 25%%: %.3f\n",
            mean(bias_med_bins > .10), mean(bias_med_bins > .25)))
cat("published vs posterior MEAN of M, per-pair median relative bias:\n"); print(round(q(P$relbias_mean_median), 4))
cat(sprintf("  per-pair fraction of bins where published < mean(M), median over pairs: %.3f\n", median(P$frac_mean_neg)))
cat("Jensen, min of posterior means vs mean of M, per-pair median relative bias:\n"); print(round(q(P$jensen_means_median), 4))
cat("pattern: Pearson cor over bins between published field and median(M), per pair:\n"); print(round(q(P$cor_pattern_med), 4))
cat("pattern: Spearman, same:\n"); print(round(q(P$spear_pattern_med), 4))
cat("pattern: Pearson cor between published field and mean(M):\n"); print(round(q(P$cor_pattern_mean), 4))
cat("where is the bias largest? Spearman cor over bins, per pair:\n")
cat("  with standardized partner gap |mu_L - mu_R| / sqrt(s_L^2 + s_R^2):\n"); print(round(q(P$cor_bias_zgap), 3))
cat("  with log-scale posterior sd of the limiting partner:\n"); print(round(q(P$cor_bias_slim), 3))
