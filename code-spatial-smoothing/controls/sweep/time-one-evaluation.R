## Time a single feature-depth evaluation (cluster + ELSA) at several depths, so the
## random-order sweep can be costed from measurement rather than guesswork.
##
## Uses the same clustering recipe as the R2.8 controls (control-A-lr-shuffle.R), which
## substitutes RcppHNSW + igraph for Seurat's neighbour backend: scale, clip at 10, no PCA,
## k = 20 KNN in the full feature space, SNN Jaccard pruned at 1/15, Louvain at resolution 1.
## Single-threaded and seeded, so runs are exactly reproducible.

suppressMessages({library(Matrix); library(RcppHNSW); library(igraph)
                  library(elsa); library(raster)})

D <- readRDS("m-est-min-matrix.rds")

## drop the four Dropbox conflicted-copy duplicates: 1481 rows -> 1477 distinct mechanisms
dup <- duplicated(D$counts)
message("duplicate mechanism rows dropped: ", sum(dup))
M   <- t(D$counts[!dup, , drop = FALSE])      # bins x mechanisms
rk  <- D$rank[!dup]
ord <- order(rk)                              # published sweep order: decreasing variance
xy  <- D$xy
message("matrix: ", nrow(M), " bins x ", ncol(M), " mechanisms")

cluster_labels <- function(M, k = 20, res = 1, seed = 42) {
  set.seed(seed)
  Xs <- scale(M); Xs[is.na(Xs)] <- 0; Xs[Xs > 10] <- 10; Xs[Xs < -10] <- -10
  idx <- hnsw_knn(Xs, k = k, distance = "euclidean", n_threads = 1)$idx
  n <- nrow(Xs)
  A <- sparseMatrix(i = rep(seq_len(n), each = k), j = as.vector(t(idx)), x = 1, dims = c(n, n))
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

mean_elsa <- function(labs, d = 50) {
  r <- rasterFromXYZ(data.frame(x = xy$x, y = xy$y, z = as.numeric(labs)))
  mean(elsa(r, d = d, categorical = TRUE)@data@values, na.rm = TRUE)
}

depths <- c(10, 150, 650, 1300, ncol(M))
res <- data.frame()
for (n in depths) {
  cols <- ord[seq_len(n)]
  t1 <- Sys.time(); labs <- cluster_labels(M[, cols, drop = FALSE]); t2 <- Sys.time()
  e   <- mean_elsa(labs);                                            t3 <- Sys.time()
  row <- data.frame(depth = n,
                    clust_s = round(as.numeric(difftime(t2, t1, units = "secs")), 1),
                    elsa_s  = round(as.numeric(difftime(t3, t2, units = "secs")), 1),
                    domains = length(unique(labs)),
                    elsa_d50 = round(e, 4))
  res <- rbind(res, row)
  message(sprintf("depth %5d | cluster %6.1f s | elsa %5.1f s | %2d domains | ELSA(d50) %.4f",
                  row$depth, row$clust_s, row$elsa_s, row$domains, row$elsa_d50))
  flush.console()
}
write.csv(res, "timing-variance-ordered.csv", row.names = FALSE)
message("\ntotal per-evaluation seconds by depth:")
print(transform(res, total_s = clust_s + elsa_s))
