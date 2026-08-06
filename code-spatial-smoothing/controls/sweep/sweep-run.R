## Random-order feature-depth sweep.
##
## Separates "how many interaction fields enter the clustering" from "which ones". The
## published sweep (Sup7, Fig 6) admits features in decreasing standardised variance, so at
## depth N it always holds the same N features and the depth effect is confounded with
## feature identity. Three arms:
##
##   ordered   decreasing variance, as published
##   reverse   increasing variance, the opposite extreme
##   random    uniform random subsets, one draw per array task
##
## Clustering follows the R2.8 controls recipe (control-A-lr-shuffle.R): scale, clip at 10,
## no PCA, k = 20 KNN in the full feature space, SNN Jaccard pruned at 1/15, Louvain at
## resolution 1. Single-threaded with a FIXED clustering seed, so between-run variation is
## attributable to the feature subset and not to clustering RNG.
##
## Runs either as a SLURM array task (one task per draw) or locally over a range:
##   SLURM:  SLURM_ARRAY_TASK_ID=7 Rscript sweep-run.R
##   local:  Rscript sweep-run.R 0 101 8          # first last ncores
##
## Each task writes one small CSV and overwrites it on restart, so --requeue is safe.

suppressMessages({
  library(Matrix); library(RcppHNSW); library(igraph); library(elsa); library(raster)
})

ROOT    <- Sys.getenv("ASOMICS_SWEEP_ROOT", unset = ".")
MATRIX  <- file.path(ROOT, "m-est-min-matrix.rds")
OUTDIR  <- file.path(ROOT, "results")
CLUSTER_SEED <- 42L
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

DESIGN_DEPTHS <- c(10, 50, 150, 350, 650, 1000, 1300)
N_RANDOM      <- 99L

## ---------------------------------------------------------------- data, once per process
D   <- readRDS(MATRIX)
dup <- duplicated(D$counts)                 # 32 rows: 4 conflicted copies + 28 operator collapses
M   <- t(D$counts[!dup, , drop = FALSE])    # bins x mechanisms
rk  <- D$rank[!dup]
xy  <- D$xy
NF  <- ncol(M)
ORD <- order(rk)                            # published order: decreasing variance
REV <- rev(ORD)
REFERENCE_DEPTHS <- unique(c(10, seq(50, NF - 50, 50), NF))

## ---------------------------------------------------------------- pipeline
cluster_labels <- function(X, k = 20, res = 1) {
  set.seed(CLUSTER_SEED)
  Xs <- scale(X); Xs[is.na(Xs)] <- 0; Xs[Xs > 10] <- 10; Xs[Xs < -10] <- -10
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

elsa_at <- function(labs, ds = c(50, 100)) {
  r <- rasterFromXYZ(data.frame(x = xy$x, y = xy$y, z = as.numeric(labs)))
  vapply(ds, function(d) mean(elsa(r, d = d, categorical = TRUE)@data@values, na.rm = TRUE), 0)
}

evaluate <- function(cols, arm, draw, depth) {
  t0 <- Sys.time()
  labs <- cluster_labels(M[, cols, drop = FALSE])
  e    <- elsa_at(labs)
  data.frame(arm = arm, draw = draw, depth = depth,
             domains = length(unique(labs)),
             elsa_d50 = e[1], elsa_d100 = e[2],
             seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
}

## ---------------------------------------------------------------- task table
## 0 ordered (design depths) | 1 reverse (design depths)
## 2 ordered (reference curve, all depths) | 3..(2+N_RANDOM) random draws
task_spec <- function(id) {
  if (id == 0L) return(list(arm = "ordered",   draw = 0L, depths = DESIGN_DEPTHS))
  if (id == 1L) return(list(arm = "reverse",   draw = 0L, depths = DESIGN_DEPTHS))
  if (id == 2L) return(list(arm = "reference", draw = 0L, depths = REFERENCE_DEPTHS))
  list(arm = "random", draw = id - 2L, depths = DESIGN_DEPTHS)
}
MAX_TASK <- 2L + N_RANDOM

run_task <- function(id) {
  s <- task_spec(id)
  out <- file.path(OUTDIR, sprintf("sweep-%s-%03d.csv", s$arm, s$draw))
  rows <- lapply(s$depths, function(n) {
    cols <- switch(s$arm,
      ordered   = ORD[seq_len(n)],
      reference = ORD[seq_len(n)],
      reverse   = REV[seq_len(n)],
      random    = { set.seed(1000L + s$draw); sample.int(NF, n) })
    r <- evaluate(cols, s$arm, s$draw, n)
    message(sprintf("[task %3d] %-9s draw %3d depth %5d -> %2d domains, ELSA(d50) %.4f, %.1f s",
                    id, s$arm, s$draw, n, r$domains, r$elsa_d50, r$seconds))
    r
  })
  res <- do.call(rbind, rows)
  write.csv(res, out, row.names = FALSE)   # idempotent: same task always rewrites this file
  invisible(res)
}

## ---------------------------------------------------------------- entry point
sid <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
if (nzchar(sid)) {
  run_task(as.integer(sid))
} else {
  a <- commandArgs(trailingOnly = TRUE)
  first <- if (length(a) >= 1) as.integer(a[1]) else 0L
  last  <- if (length(a) >= 2) as.integer(a[2]) else MAX_TASK
  ncore <- if (length(a) >= 3) as.integer(a[3]) else 1L
  ids <- seq(first, min(last, MAX_TASK))
  message("running tasks ", first, "..", max(ids), " on ", ncore, " core(s); ",
          NF, " distinct mechanisms")
  if (ncore > 1L) {
    parallel::mclapply(ids, run_task, mc.cores = ncore, mc.preschedule = FALSE)
  } else {
    lapply(ids, run_task)
  }
  message("done")
}
