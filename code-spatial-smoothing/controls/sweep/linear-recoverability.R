## How much of the min field is recoverable by a LINEAR function of its two parents?
##
## This is the direct pre-check for question (b): the published pipeline builds a
## Euclidean kNN graph on z-scored features with no PCA, so its metric is linear in
## whatever features it is handed. If min(L,R) is well approximated by a linear
## combination of L and R, a pipeline given L and R separately can already express
## almost everything the min field carries, and the explicit gate buys little.
## If the linear fit is poor, the gate encodes something the metric cannot build.
##
## Per mechanism: R^2 of lm(min ~ L + R), compared with r^2 of the dominant parent alone.

mn <- readRDS("m-est-min-matrix.rds"); gf <- readRDS("gene-field-matrix.rds")
e <- new.env(); load("fantom-hierarchy.Robj", envir = e); fan <- get("output", e)
M <- mn$counts; G <- gf$G
fan$key <- gsub("_", "-", fan$MECHANISM)
idx <- match(rownames(M), fan$key); M <- M[!is.na(idx), ]; idx <- idx[!is.na(idx)]
gup <- toupper(colnames(G))
li <- match(fan$LIGAND[idx], gup); ri <- match(fan$RECEPTOR[idx], gup)
ok <- !is.na(li) & !is.na(ri); M <- M[ok, ]; li <- li[ok]; ri <- ri[ok]
n <- nrow(M)

R2_both <- numeric(n); r2_dom <- numeric(n); dominance <- numeric(n)
for (i in seq_len(n)) {
  L <- G[, li[i]]; R <- G[, ri[i]]; m <- M[i, ]
  fl <- mean(L <= R); dominance[i] <- max(fl, 1 - fl)
  D <- if (fl >= 0.5) L else R
  r2_dom[i] <- suppressWarnings(cor(m, D))^2
  fit <- .lm.fit(cbind(1, L, R), m)
  R2_both[i] <- 1 - sum(fit$residuals^2) / sum((m - mean(m))^2)
}
gain <- R2_both - r2_dom

f <- function(x) format(round(quantile(x, c(.1,.25,.5,.75,.9), na.rm=TRUE), 4), nsmall=4)
cat("R^2 of the BEST LINEAR fit of min on (L, R)\n")
cat("  10/25/50/75/90 pct: ", paste(f(R2_both), collapse="  "), "\n")
cat("  median: ", round(median(R2_both),4), "\n", sep="")
for (t in c(.90,.95,.99)) cat("  R^2 > ", t, ": ", sum(R2_both > t), " of ", n,
  " (", round(100*mean(R2_both > t),1), "%)\n", sep="")
cat("\n  UNEXPLAINED by any linear use of both parents (1 - R^2), median: ",
    round(median(1 - R2_both), 4), "\n", sep="")
cat("  mechanisms where linear leaves >10% of variance unexplained: ",
    sum(R2_both < .90), " (", round(100*mean(R2_both < .90),1), "%)\n\n", sep="")

cat("r^2 of the DOMINANT PARENT alone, median: ", round(median(r2_dom),4), "\n", sep="")
cat("gain from adding the second parent linearly, median: ", round(median(gain),4), "\n\n", sep="")

hi <- dominance > .95
cat("split by dominance:\n")
cat("  dominance > 0.95 (n=", sum(hi), "): median linear R^2 = ",
    round(median(R2_both[hi]),4), "\n", sep="")
cat("  dominance <= 0.95 (n=", sum(!hi), "): median linear R^2 = ",
    round(median(R2_both[!hi]),4), "\n", sep="")

write.csv(data.frame(mechanism=rownames(M), dominance=dominance,
  r2_dominant_parent=r2_dom, R2_linear_both=R2_both, gain=gain),
  "min-linear-recoverability.csv", row.names=FALSE)
cat("\nwrote min-linear-recoverability.csv\n")
