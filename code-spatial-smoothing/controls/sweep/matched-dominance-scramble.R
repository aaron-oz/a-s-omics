## Matched-dominance scramble: is the minimum's nonlinear content pairing-specific?
##
## min(L,R) carries structure that no linear function of L and R recovers (median
## R^2 = 0.555, measured in linear-recoverability.R). The amount of that structure is
## driven mainly by DOMINANCE, the fraction of bins at which the same partner supplies
## the minimum: switching creates the nonlinearity, and a mechanism with dominance 1.0
## has none at all (min is exactly one parent).
##
## So the fair null is not a random partner, it is a random partner MATCHED ON
## DOMINANCE. That holds the main driver constant and isolates the pairing.
##
## Two statistics per mechanism, computed for the true pairing and its matched scramble:
##   linear R^2       R^2 of lm(min ~ L + R). 1 - R^2 is the nonlinear content.
##   switch coherence fraction of 4-neighbour adjacent bin pairs whose limiting partner
##                    is the same. High = switching follows smooth spatial domains;
##                    low = the operator flips partner bin-to-bin, i.e. the extra
##                    structure is high-frequency artifact rather than anatomy.

set.seed(11)
NCAND <- 30

mn <- readRDS("m-est-min-matrix.rds"); gf <- readRDS("gene-field-matrix.rds")
e <- new.env(); load("fantom-hierarchy.Robj", envir = e); fan <- get("output", e)
M <- mn$counts; G <- gf$G; XY <- gf$xy

fan$key <- gsub("_", "-", fan$MECHANISM)
idx <- match(rownames(M), fan$key); M <- M[!is.na(idx), ]; idx <- idx[!is.na(idx)]
gup <- toupper(colnames(G))
li <- match(fan$LIGAND[idx], gup); ri <- match(fan$RECEPTOR[idx], gup)
ok <- !is.na(li) & !is.na(ri); M <- M[ok, ]; li <- li[ok]; ri <- ri[ok]
n <- nrow(M); cat("mechanisms: ", n, "\n", sep = "")

## exclusion set: every real FANTOM5 pair, in gene-index space, both directions
alL <- match(fan$LIGAND, gup); alR <- match(fan$RECEPTOR, gup)
keep <- !is.na(alL) & !is.na(alR)
REAL <- unique(c(paste(alL[keep], alR[keep], sep = "|"),
                 paste(alR[keep], alL[keep], sep = "|")))
cat("real FANTOM5 pairs excluded (both directions): ", length(REAL), "\n", sep = "")

## spatial grid for the switch-coherence statistic
ux <- sort(unique(XY[, 1])); uy <- sort(unique(XY[, 2]))
gx <- match(XY[, 1], ux);    gy <- match(XY[, 2], uy)
cat("grid: ", length(ux), " x ", length(uy), "\n\n", sep = "")
cell <- matrix(NA_integer_, length(ux), length(uy)); cell[cbind(gx, gy)] <- seq_len(nrow(XY))

switch_coherence <- function(mask) {         # mask TRUE where the ligand is limiting
  A <- matrix(NA, length(ux), length(uy)); A[cbind(gx, gy)] <- mask
  h <- A[-1, , drop = FALSE] == A[-nrow(A), , drop = FALSE]
  v <- A[, -1, drop = FALSE] == A[, -ncol(A), drop = FALSE]
  mean(c(h, v), na.rm = TRUE)
}

stats_for <- function(L, R) {
  m  <- pmin(L, R)
  fl <- mean(L <= R)
  if (sd(m) == 0) return(c(dom = max(fl, 1-fl), R2 = NA, sc = NA))
  fit <- .lm.fit(cbind(1, L, R), m)
  c(dom = max(fl, 1 - fl),
    R2  = 1 - sum(fit$residuals^2) / sum((m - mean(m))^2),
    sc  = switch_coherence(L <= R))
}

rpool <- sort(unique(ri))                    # candidate partners: the fitted receptors
TRUE_ <- matrix(NA_real_, n, 3); NULL_ <- matrix(NA_real_, n, 3)
matched_gene <- integer(n); dom_gap <- numeric(n)

t0 <- Sys.time()
for (i in seq_len(n)) {
  L <- G[, li[i]]; R <- G[, ri[i]]
  TRUE_[i, ] <- stats_for(L, R)
  d_true <- TRUE_[i, 1]

  cand <- rpool[rpool != li[i] & rpool != ri[i]]
  cand <- cand[!(paste(li[i], cand, sep = "|") %in% REAL)]
  cand <- if (length(cand) > NCAND) sample(cand, NCAND) else cand
  Rc <- G[, cand, drop = FALSE]
  fl <- colMeans(L <= Rc); dcand <- pmax(fl, 1 - fl)
  j <- which.min(abs(dcand - d_true))
  matched_gene[i] <- cand[j]; dom_gap[i] <- abs(dcand[j] - d_true)
  NULL_[i, ] <- stats_for(L, G[, cand[j]])
  if (i %% 300 == 0) cat("  ", i, "/", n, "\n", sep = "")
}
cat("elapsed: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min\n\n", sep = "")

sm <- function(lab, x) cat(sprintf("  %-22s median %.4f   mean %.4f   IQR %.4f-%.4f\n",
  lab, median(x, na.rm=TRUE), mean(x, na.rm=TRUE),
  quantile(x,.25,na.rm=TRUE), quantile(x,.75,na.rm=TRUE)))

cat("MATCHING QUALITY (|dominance_null - dominance_true|)\n")
cat("  median gap ", round(median(dom_gap), 5), ", 90th pct ",
    round(quantile(dom_gap, .9), 5), ", max ", round(max(dom_gap), 5), "\n", sep="")
cat("  dominance, true: median ", round(median(TRUE_[,1]),4),
    " | null: median ", round(median(NULL_[,1]),4), "\n\n", sep="")

cat("LINEAR R^2 of min ~ L + R  (lower = more nonlinear content)\n")
sm("true pairing", TRUE_[,2]); sm("matched scramble", NULL_[,2])
d <- TRUE_[,2] - NULL_[,2]
cat(sprintf("  paired difference: median %.4f  mean %.4f  sd %.4f\n", median(d,na.rm=TRUE), mean(d,na.rm=TRUE), sd(d,na.rm=TRUE)))
cat(sprintf("  Wilcoxon signed-rank p = %.3g ; true more nonlinear in %.1f%% of mechanisms\n\n",
    wilcox.test(TRUE_[,2], NULL_[,2], paired=TRUE)$p.value, 100*mean(d < 0, na.rm=TRUE)))

cat("SWITCH COHERENCE (fraction of adjacent bins with the same limiting partner)\n")
sm("true pairing", TRUE_[,3]); sm("matched scramble", NULL_[,3])
d2 <- TRUE_[,3] - NULL_[,3]
cat(sprintf("  paired difference: median %.4f  mean %.4f  sd %.4f\n", median(d2,na.rm=TRUE), mean(d2,na.rm=TRUE), sd(d2,na.rm=TRUE)))
cat(sprintf("  Wilcoxon signed-rank p = %.3g ; true more coherent in %.1f%% of mechanisms\n",
    wilcox.test(TRUE_[,3], NULL_[,3], paired=TRUE)$p.value, 100*mean(d2 > 0, na.rm=TRUE)))
cat(sprintf("  reference: a spatially random mask of the same dominance would sit near %.4f\n\n",
    mean(TRUE_[,1]^2 + (1-TRUE_[,1])^2)))

write.csv(data.frame(mechanism=rownames(M), ligand=colnames(G)[li], receptor=colnames(G)[ri],
  matched_partner=colnames(G)[matched_gene], dom_gap=dom_gap,
  dom_true=TRUE_[,1], R2_true=TRUE_[,2], switchcoh_true=TRUE_[,3],
  dom_null=NULL_[,1], R2_null=NULL_[,2], switchcoh_null=NULL_[,3]),
  "matched-dominance-scramble.csv", row.names=FALSE)
cat("wrote matched-dominance-scramble.csv\n")
