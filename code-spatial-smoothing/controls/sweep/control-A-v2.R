## Control A, second attempt: ligand-receptor PAIR-SHUFFLE null, done properly.
##
## The original Control A was withdrawn for two independent reasons:
##   1. it compared single runs of a clustering step that did not reproduce between runs,
##      so the effect was smaller than the pipeline's own noise; and
##   2. its scramble was naive, so a receptor permutation recreated a mean of 73.5 genuine
##      FANTOM5 mechanisms out of 1477 (5.0%), meaning the "null" was ~5% real pairs.
##
## Both are now fixable. The clustering is exactly reproducible at one thread with a fixed
## seed, so single runs are interpretable. And the scramble below is exclusion-aware: after
## permuting receptors it repairs every pair that coincides with a real FANTOM5 mechanism,
## and every self-pair, by targeted swaps, so the null contains exactly zero real pairs.
##
## Question: does the TRUE ligand-receptor pairing yield more spatially coherent tissue
## domains than scrambled pairings built from the identical set of gene fields, with the
## same dimensionality and the same clustering? If not, the interaction fields are recovering
## tissue architecture from the gene fields rather than from the pairing.
##
## One array task = one scramble. Task 0 = the true pairing.
##   SLURM:  SLURM_ARRAY_TASK_ID=0 Rscript control-A-v2.R
##   local:  Rscript control-A-v2.R 0 5        # first last

suppressMessages({
  library(data.table); library(RcppHNSW); library(igraph); library(Matrix)
  library(elsa); library(raster)
})

ROOT   <- Sys.getenv("ASOMICS_CTRLA_ROOT", unset = ".")
OUTDIR <- file.path(ROOT, "results-ctrlA")
CLUSTER_SEED <- 42L
N_SCRAMBLE   <- 99L
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------- inputs
GG <- readRDS(file.path(ROOT, "gene-field-matrix.rds"))
G  <- GG$G; xy <- GG$xy                     # locations x genes, posterior medians
GENE <- toupper(colnames(G))

e <- new.env(); load(file.path(ROOT, "fantom-hierarchy.Robj"), envir = e)
hier <- as.data.table(get(ls(e)[1], envir = e))
mech <- unique(hier[, .(L = LIGAND, R = RECEPTOR)])
mech <- mech[L %in% GENE & R %in% GENE]
mech[, li := match(L, GENE)][, ri := match(R, GENE)]
mech <- mech[!is.na(li) & !is.na(ri)]

## every real FANTOM5 pair, including ones whose partners are not both fitted, so the
## exclusion is against the full database rather than only the modelled subset
REAL <- unique(paste(hier$LIGAND, hier$RECEPTOR, sep = "|"))
message(sprintf("mechanisms with both partners fitted: %d; real pairs in FANTOM5: %d",
                nrow(mech), length(REAL)))

## ---------------------------------------------------------------- pipeline
cluster_labels <- function(M, k = 20, res = 1) {
  set.seed(CLUSTER_SEED)
  Xs <- scale(M); Xs[is.na(Xs)] <- 0; Xs[Xs > 10] <- 10; Xs[Xs < -10] <- -10
  idx <- hnsw_knn(Xs, k = k, distance = "euclidean", n_threads = 1)$idx
  n <- nrow(Xs)
  A  <- sparseMatrix(i = rep(seq_len(n), each = k), j = as.vector(t(idx)), x = 1, dims = c(n, n))
  ov <- as(A %*% t(A), "TsparseMatrix")
  jac <- ov@x / (2 * k - ov@x); keep <- jac > (1/15) & ov@i != ov@j
  g <- simplify(graph_from_data_frame(
        data.frame(from = ov@i[keep] + 1, to = ov@j[keep] + 1, weight = jac[keep]),
        directed = FALSE))
  memb <- rep(NA_integer_, n); cl <- cluster_louvain(g, resolution = res)
  memb[as.integer(V(g)$name)] <- membership(cl)
  if (any(is.na(memb))) memb[is.na(memb)] <- max(memb, na.rm = TRUE) + seq_len(sum(is.na(memb)))
  memb
}

elsa_mean <- function(labels, dd) {
  r <- rasterFromXYZ(cbind(xy[, 1], xy[, 2], as.integer(labels)))
  mean(elsa(r, d = dd, categorical = TRUE)@data@values, na.rm = TRUE)
}

## Permute receptors, then repair coincidences by targeted swaps. A pair is bad if it is a
## real FANTOM5 mechanism or a self-pair. Preserves the receptor multiset exactly, so the
## marginal field statistics of the null are identical to the observed.
scramble <- function(li, ri, seed) {
  set.seed(seed)
  n <- length(li)
  ri <- sample(ri)
  is_bad <- function(a, b) (paste(GENE[a], GENE[b], sep = "|") %in% REAL) | (a == b)
  for (sweep in 1:200) {
    bad <- which(is_bad(li, ri))
    if (!length(bad)) break
    for (i in bad) {
      for (attempt in 1:100) {
        j <- sample.int(n, 1L)
        if (j == i) next
        if (!is_bad(li[i], ri[j]) && !is_bad(li[j], ri[i])) {
          tmp <- ri[i]; ri[i] <- ri[j]; ri[j] <- tmp; break
        }
      }
    }
  }
  list(ri = ri, remaining_real = sum(is_bad(li, ri)))
}

evaluate <- function(li, ri, tag, remaining = 0L) {
  t0 <- Sys.time()
  M  <- pmin(G[, li, drop = FALSE], G[, ri, drop = FALSE])   # common-minimum convolution
  ndist <- sum(!duplicated(t(M)))                            # distinct fields after collapse
  labs <- cluster_labels(M)
  data.frame(tag = tag, n_mech = length(li), n_distinct = ndist,
             remaining_real_pairs = remaining,
             n_clusters = length(unique(labs)),
             elsa_d50 = elsa_mean(labs, 50), elsa_d100 = elsa_mean(labs, 100),
             seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
}

run_task <- function(id) {
  if (id == 0L) {
    res <- evaluate(mech$li, mech$ri, "observed", 0L)
    out <- file.path(OUTDIR, "ctrlA-observed.csv")
  } else {
    s <- scramble(mech$li, mech$ri, seed = 5000L + id)
    res <- evaluate(mech$li, s$ri, sprintf("scramble%03d", id), s$remaining_real)
    out <- file.path(OUTDIR, sprintf("ctrlA-scramble-%03d.csv", id))
  }
  write.csv(res, out, row.names = FALSE)
  message(sprintf("[%s] %d mech, %d distinct, %d real pairs left, %d clusters, ELSA(d50) %.4f, %.0f s",
                  res$tag, res$n_mech, res$n_distinct, res$remaining_real_pairs,
                  res$n_clusters, res$elsa_d50, res$seconds))
  invisible(res)
}

sid <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
if (nzchar(sid)) {
  run_task(as.integer(sid))
} else {
  a <- commandArgs(trailingOnly = TRUE)
  first <- if (length(a) >= 1) as.integer(a[1]) else 0L
  last  <- if (length(a) >= 2) as.integer(a[2]) else N_SCRAMBLE
  for (i in seq(first, last)) run_task(i)
}
