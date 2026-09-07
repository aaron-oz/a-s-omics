## Does the common minimum return a re-scaled copy of one parent gene field?
##
## Motivation: the min is taken on raw counts, so it is decided by which partner is
## globally dimmer. Seurat then z-scores every feature before building the Euclidean
## kNN graph, which deletes magnitude. If a min field is highly correlated with one
## parent, then after scaling it is a near-duplicate of a column the 912-gene input
## set already contains, and it contributes almost nothing new to the metric.
##
## Pearson r is invariant to the z-scoring, so r computed here IS the post-scaling
## similarity. For unit-variance vectors ||a-b||^2 = 2n(1-r), so r is sufficient.

t0 <- Sys.time()
mn <- readRDS("m-est-min-matrix.rds")
gf <- readRDS("gene-field-matrix.rds")
e <- new.env(); load("fantom-hierarchy.Robj", envir = e); fan <- get("output", e)

M <- mn$counts                      # 1481 mechanisms x 26000 bins
G <- gf$G                           # 26000 bins x 912 genes

## --- align bins -------------------------------------------------------------
key_g <- paste(gf$xy[, 1], gf$xy[, 2], sep = ".")
stopifnot(identical(key_g, colnames(M)))
cat("bin alignment: exact, ", ncol(M), " bins\n", sep = "")

## --- drop the 4 duplicated mechanism rows (Dropbox conflicted copies) --------
dup <- duplicated(rownames(M))
cat("duplicated mechanism names dropped: ", sum(dup), "\n", sep = "")
M <- M[!dup, , drop = FALSE]

## --- map mechanism -> ligand, receptor --------------------------------------
fan$key <- gsub("_", "-", fan$MECHANISM)
idx <- match(rownames(M), fan$key)
cat("mechanisms matched to FANTOM: ", sum(!is.na(idx)), " of ", nrow(M), "\n", sep = "")
M <- M[!is.na(idx), , drop = FALSE]; idx <- idx[!is.na(idx)]

gup <- toupper(colnames(G))
li <- match(fan$LIGAND[idx],   gup)
ri <- match(fan$RECEPTOR[idx], gup)
ok <- !is.na(li) & !is.na(ri)
cat("mechanisms with both partners fitted: ", sum(ok), "\n", sep = "")
M <- M[ok, , drop = FALSE]; li <- li[ok]; ri <- ri[ok]
n <- nrow(M); cat("analysis set: ", n, " mechanisms\n\n", sep = "")

## --- validate: is the stored field really pmin(L, R)? -----------------------
set.seed(1); chk <- sample.int(n, 25)
mx <- max(sapply(chk, function(i)
  max(abs(M[i, ] - pmin(G[, li[i]], G[, ri[i]])))))
cat("VALIDATION max |stored - pmin(L,R)| over 25 mechanisms: ",
    format(mx, digits = 3), "\n\n", sep = "")

## --- dominance and parent correlation ---------------------------------------
zs <- function(v) { s <- sd(v); if (s == 0) rep(0, length(v)) else (v - mean(v)) / s }

frac_lig <- numeric(n); r_dom <- numeric(n); r_oth <- numeric(n); dom <- character(n)
for (i in seq_len(n)) {
  L <- G[, li[i]]; R <- G[, ri[i]]; m <- M[i, ]
  fl <- mean(L <= R); frac_lig[i] <- fl
  rl <- suppressWarnings(cor(m, L)); rr <- suppressWarnings(cor(m, R))
  if (fl >= 0.5) { dom[i] <- "ligand"; r_dom[i] <- rl; r_oth[i] <- rr }
  else           { dom[i] <- "receptor"; r_dom[i] <- rr; r_oth[i] <- rl }
}
dominance <- pmax(frac_lig, 1 - frac_lig)

q <- function(x) format(round(quantile(x, c(.1,.25,.5,.75,.9), na.rm=TRUE), 4), nsmall=4)
cat("DOMINANCE (fraction of bins where the same partner supplies the minimum)\n")
cat("  deciles/quartiles/median: ", paste(q(dominance), collapse = "  "), "\n")
cat("  >0.95 at: ", sum(dominance > .95), " (", round(100*mean(dominance > .95),1), "%)\n", sep="")
cat("  >0.90 at: ", sum(dominance > .90), " (", round(100*mean(dominance > .90),1), "%)\n", sep="")
cat("  ==1.00  : ", sum(dominance == 1), " (min field is exactly one parent field)\n\n", sep="")

cat("CORRELATION of the min field with its DOMINANT parent\n")
cat("  quantiles: ", paste(q(r_dom), collapse = "  "), "\n")
cat("  median overall: ", round(median(r_dom, na.rm=TRUE), 4), "\n", sep="")
for (thr in c(.90, .95, .99)) cat("  r > ", thr, ": ", sum(r_dom > thr, na.rm=TRUE),
    " (", round(100*mean(r_dom > thr, na.rm=TRUE), 1), "%)\n", sep="")
cat("\n  restricted to the ", sum(dominance > .95), " mechanisms with dominance > 0.95:\n", sep="")
s <- dominance > .95
cat("    median r with dominant parent: ", round(median(r_dom[s], na.rm=TRUE), 4), "\n", sep="")
for (thr in c(.95, .99)) cat("    r > ", thr, ": ", sum(r_dom[s] > thr, na.rm=TRUE),
    " of ", sum(s), "\n", sep="")
cat("\n  correlation with the OTHER (non-dominant) parent, median: ",
    round(median(r_oth, na.rm=TRUE), 4), "\n\n", sep="")

## --- strongest version: nearest neighbour among ALL 912 input gene fields ----
cat("NEAREST INPUT GENE FIELD (max |r| against all 912 L/R gene fields)\n")
Mz <- t(apply(M, 1, zs))            # n x 26000
Gz <- apply(G, 2, zs)               # 26000 x 912
C  <- (Mz %*% Gz) / (ncol(M) - 1)   # n x 912 Pearson
best <- apply(abs(C), 1, max)
bestg <- colnames(G)[apply(abs(C), 1, which.max)]
cat("  quantiles: ", paste(q(best), collapse = "  "), "\n")
cat("  median: ", round(median(best), 4), "\n", sep="")
for (thr in c(.90, .95, .99)) cat("  max|r| > ", thr, ": ", sum(best > thr),
    " of ", n, " (", round(100*mean(best > thr), 1), "%)\n", sep="")
cat("  is the nearest gene one of its own two partners? ",
    sum(bestg == colnames(G)[li] | bestg == colnames(G)[ri]), " of ", n, "\n\n", sep="")

out <- data.frame(mechanism = rownames(M), ligand = colnames(G)[li],
                  receptor = colnames(G)[ri], frac_from_ligand = frac_lig,
                  dominance = dominance, dominant = dom, r_dominant = r_dom,
                  r_other = r_oth, max_r_any_input_gene = best, nearest_gene = bestg)
write.csv(out, "min-vs-parent-correlation.csv", row.names = FALSE)
cat("wrote min-vs-parent-correlation.csv (", nrow(out), " rows)\n", sep="")
cat("elapsed: ", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), " s\n", sep="")
