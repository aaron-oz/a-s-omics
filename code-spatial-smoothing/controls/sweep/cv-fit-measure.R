## Cross-validated predictive fit: the candidate primary measure, plus its calibration.
##
## Scores a partition by how well it predicts held-out gene expression in bins it was not
## fitted on. A spurious split gives two pieces with the same true profile, so it adds a
## parameter without improving out-of-sample prediction. A real split gives two pieces with
## different profiles, so prediction improves. Cross-validation does the penalising, which
## removes the sibling definition, the spatial null, the alpha, the multiplicity correction
## and the minimum domain size that the previous candidate measure needed.
##
## Model. Bin j in domain d, gene g:  y_gj ~ Poisson(N_j * mu_gd), with N_j the bin's total
## count. This is the same offset structure the paper's own model uses. Domain rates are
## shrunk toward the global rate by one bin's worth of depth, which also handles domains with
## no training bins in a fold.
##
## Folds are contiguous spatial blocks, not scattered bins: with scattered bins a test bin is
## surrounded by training bins and spatial smoothness predicts it regardless of the partition.
## Block side is 20 bins, comfortably above the fitted correlation length (the model's range
## prior puts 95% of its mass below 500 coordinate units, which is 10 bins at bin-50).

suppressMessages({library(Matrix); library(data.table)})
set.seed(1)

ROOT <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
HERE <- file.path(ROOT, "code-spatial-smoothing/controls/sweep")

H  <- readRDS(file.path(HERE, "heldout-gene-matrix.rds"))
Y  <- H$counts                                   # genes x bins, sparse
xy <- do.call(rbind, strsplit(H$bins, "\\."))
xy <- data.table(x = as.numeric(xy[, 1]), y = as.numeric(xy[, 2]))

## the pipeline's own gene filter, applied to the held-out set
det <- Matrix::rowSums(Y > 0); two <- Matrix::rowSums(Y >= 2)
keep <- det > 0.0025 * ncol(Y) & two > 0   # the Sup1 filter, computed sparsely
Y <- Y[keep, ]
message(sprintf("held-out genes after the pipeline's own filter: %d of %d", nrow(Y), length(keep)))

N <- Matrix::colSums(Y)                          # per-bin depth over held-out genes
BLOCK <- 20; NFOLD <- 10
blk <- paste(xy$x %/% (BLOCK * 50), xy$y %/% (BLOCK * 50))
ub  <- unique(blk); fold_of_block <- setNames(sample(rep_len(1:NFOLD, length(ub))), ub)
fold <- fold_of_block[blk]
message(sprintf("%d spatial blocks of %dx%d bins across %d folds", length(ub), BLOCK, BLOCK, NFOLD))

glob <- Matrix::rowSums(Y) / sum(N)              # global per-gene rate
SHRINK <- median(N)                              # one bin's worth of depth

## Poisson deviance of the test bins under a set of domain rates, computed from the nonzeros
## plus row and column sums rather than by materialising genes x bins.
dev_test <- function(Yte, Nte, dom_te, MU) {
  S <- Matrix::colSums(MU)                       # total rate per domain
  sum_lambda <- sum(Nte * S[dom_te])
  tri <- as(Yte, "TsparseMatrix")
  gi <- tri@i + 1L; bj <- tri@j + 1L; yv <- tri@x
  lam <- Nte[bj] * MU[cbind(gi, dom_te[bj])]
  2 * (sum(yv * log(yv / lam)) - sum(yv) + sum_lambda)
}

cv_deviance_explained <- function(labels) {
  labels <- as.integer(as.factor(labels)); K <- max(labels)
  dm <- dn <- 0
  for (f in 1:NFOLD) {
    tr <- which(fold != f); te <- which(fold == f)
    Ytr <- Y[, tr, drop = FALSE]; Ntr <- N[tr]
    ind <- sparseMatrix(i = seq_along(tr), j = labels[tr], x = 1, dims = c(length(tr), K))
    Ycnt <- as.matrix(Ytr %*% ind)               # genes x domains, training counts
    Ndom <- as.numeric(Matrix::crossprod(ind, Ntr))
    MU <- sweep(Ycnt + SHRINK * glob, 2, Ndom + SHRINK, "/")
    MU0 <- matrix(glob, nrow = nrow(Y), ncol = K)
    dm <- dm + dev_test(Y[, te, drop = FALSE], N[te], labels[te], MU)
    dn <- dn + dev_test(Y[, te, drop = FALSE], N[te], labels[te], MU0)
  }
  1 - dm / dn
}

## ---- the partition the paper published --------------------------------------------------
e <- new.env()
load(file.path(ROOT, "sam_sandbox/clustering vs features/m.est.min",
                "clustering vs features experiment",
                "meta.data.feature.depth.experiment.2025-07-12.Robj"), envir = e)
md <- get(ls(e)[1], envir = e)
cc <- grep("featureClusters$", colnames(md), value = TRUE)
base <- as.integer(as.factor(md[[cc[which(as.integer(sub("featureClusters$","",cc)) == 1481)]]]))
stopifnot(length(base) == ncol(Y))
message(sprintf("published partition: %d domains over %d bins", max(base), length(base)))

ann <- as.data.frame(fread(file.path(ROOT,
        "data-outputs/mouse-embryo-raw/controls/e16-fov-annotation.csv")))
ord <- match(paste(xy$x, xy$y, sep = "."), paste(ann$x, ann$y, sep = "."))
stopifnot(!anyNA(ord))
ann <- ann[ord, ]

purity_of <- function(lab, d) {
  a <- ann$annotation[lab == d]; max(table(a)) / length(a)
}
pur <- sapply(sort(unique(base)), function(d) purity_of(base, d))
sz  <- as.integer(table(base))

## ---- surgery operations ------------------------------------------------------------------
## arbitrary split: cut a domain in two along its longer spatial axis, a boundary with no
## reason to exist
split_arbitrary <- function(lab, d) {
  i <- which(lab == d)
  v <- if (diff(range(xy$x[i])) >= diff(range(xy$y[i]))) xy$x[i] else xy$y[i]
  lab[i[v > median(v)]] <- max(lab) + 1L
  lab
}
## informative split: cut a domain along the reference annotation's own boundary inside it,
## which is a division we have external reason to believe is real
split_informative <- function(lab, d) {
  i <- which(lab == d); a <- ann$annotation[i]
  top <- names(sort(table(a), decreasing = TRUE))[1]
  if (sum(a != top) < 25) return(NULL)
  lab[i[a != top]] <- max(lab) + 1L
  lab
}
merge_two <- function(lab, d1, d2) { lab[lab == d2] <- d1; lab }
jitter_boundary <- function(lab, frac = 0.10) {
  key <- paste(xy$x, xy$y); idx <- setNames(seq_along(key), key)
  nb <- cbind(idx[paste(xy$x + 50, xy$y)], idx[paste(xy$x - 50, xy$y)],
              idx[paste(xy$x, xy$y + 50)], idx[paste(xy$x, xy$y - 50)])
  onb <- which(apply(nb, 1, function(r) { r <- r[!is.na(r)]; any(lab[r] != lab[r[1]]) |
                                           any(lab[r] != lab[which(nb[1,1] == nb[1,1])][1]) }))
  onb <- which(sapply(seq_len(nrow(nb)), function(k) {
    r <- nb[k, ]; r <- r[!is.na(r)]; any(lab[r] != lab[k]) }))
  pick <- sample(onb, floor(frac * length(onb)))
  for (k in pick) { r <- nb[k, ]; r <- r[!is.na(r)]; r <- r[lab[r] != lab[k]]
    if (length(r)) lab[k] <- lab[r[1]] }
  lab
}

## adjacency, for choosing a merge pair
adj_pairs <- function(lab) {
  key <- paste(xy$x, xy$y); idx <- setNames(seq_along(key), key)
  r <- idx[paste(xy$x + 50, xy$y)]; d <- idx[paste(xy$x, xy$y + 50)]
  p <- rbind(cbind(lab, lab[r]), cbind(lab, lab[d]))
  p <- p[complete.cases(p) & p[,1] != p[,2], , drop = FALSE]
  t(apply(p, 1, sort))
}

message("\nscoring baseline ...")
t0 <- Sys.time(); b <- cv_deviance_explained(base)
message(sprintf("  baseline deviance explained %.5f  (%.0f s)", b,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))

res <- data.frame(test = "baseline", detail = sprintf("%d domains", max(base)),
                  score = b, delta = 0, stringsAsFactors = FALSE)
add <- function(test, detail, s) res <<- rbind(res, data.frame(test = test, detail = detail,
                                    score = s, delta = s - b, stringsAsFactors = FALSE))

## negative control: arbitrary splits of the three purest domains
for (d in order(-pur)[1:3]) {
  s <- cv_deviance_explained(split_arbitrary(base, d))
  add("arbitrary split", sprintf("domain %d (n=%d, purity %.2f)", d, sz[d], pur[d]), s)
  message(sprintf("  arbitrary split of domain %2d: %+.5f", d, s - b))
}
## positive control: split the three least pure domains along the reference boundary
for (d in order(pur)[1:3]) {
  L <- split_informative(base, d); if (is.null(L)) next
  s <- cv_deviance_explained(L)
  add("informative split", sprintf("domain %d (n=%d, purity %.2f)", d, sz[d], pur[d]), s)
  message(sprintf("  informative split of domain %2d: %+.5f", d, s - b))
}
## merges of adjacent pairs
ap <- adj_pairs(base); tab <- table(paste(ap[,1], ap[,2]))
for (k in names(sort(tab, decreasing = TRUE))[1:3]) {
  pr <- as.integer(strsplit(k, " ")[[1]])
  s <- cv_deviance_explained(merge_two(base, pr[1], pr[2]))
  add("merge", sprintf("domains %d + %d", pr[1], pr[2]), s)
  message(sprintf("  merge %2d + %2d: %+.5f", pr[1], pr[2], s - b))
}
## boundary jitter
for (fr in c(0.10, 0.25)) {
  s <- cv_deviance_explained(jitter_boundary(base, fr))
  add("boundary jitter", sprintf("%.0f%% of boundary bins", 100 * fr), s)
  message(sprintf("  jitter %.0f%%: %+.5f", 100 * fr, s - b))
}

write.csv(res, file.path(HERE, "cv-fit-calibration.csv"), row.names = FALSE)
message("\n==== CALIBRATION ====")
print(transform(res, score = round(score, 5), delta = round(delta, 5)), row.names = FALSE)
