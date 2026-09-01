## Two blocking checks before the representation-comparison experiment (design v3).
##
## 1. DETERMINISM. Seurat's neighbour backend (annoy) is approximate and order-dependent. If
##    the pipeline is not reproducible given identical input and seed, then replicate spread
##    in the planned experiment conflates representation instability with clustering noise,
##    and the stability instrument means nothing.
##
## 2. NEIGHBOUR-GRAPH REUSE. FindNeighbors depends only on the feature matrix, not the
##    resolution, so sweeping FindClusters on a built graph may be nearly free. If so, the
##    same budget buys 20 replicates rather than 5, which decides whether the design can
##    support its own decision rule.
##
## Run at depth 350 rather than full depth: determinism is a property of the algorithm, and
## the reuse ratio measured here understates the saving at full depth (graph cost grows with
## features, clustering cost does not), so it is a conservative estimate. One task, ~1 h.

local({
  ul <- Sys.getenv("R_LIBS_USER")
  if (nzchar(ul) && dir.exists(ul)) .libPaths(c(ul, .libPaths()))
})
options(future.globals.maxSize = 10.0 * 1e9)
suppressMessages({library(Seurat); library(Matrix); library(mclust)})

ROOT  <- Sys.getenv("ASOMICS_SWEEP_ROOT", unset = ".")
DEPTH <- as.integer(Sys.getenv("ASOMICS_CHECK_DEPTH", "350"))
SEED  <- 42L
RES   <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0)   # the v3 resolution grid
OUT   <- file.path(ROOT, "results-checks"); dir.create(OUT, showWarnings = FALSE)

message("Seurat ", packageVersion("Seurat"), ", ggplot2 ", packageVersion("ggplot2"))

D   <- readRDS(file.path(ROOT, "m-est-min-matrix.rds"))
dup <- duplicated(D$counts)
cnt <- D$counts[!dup, , drop = FALSE]
ord <- order(D$rank[!dup])
rownames(cnt) <- make.unique(rownames(cnt)); colnames(cnt) <- paste0("bin", seq_len(ncol(cnt)))
foi <- rownames(cnt)[ord[seq_len(DEPTH)]]
message(sprintf("%d distinct mechanisms, using top %d by variance", nrow(cnt), DEPTH))

tm <- function(expr) { t0 <- Sys.time(); v <- force(expr)
  list(v = v, s = as.numeric(difftime(Sys.time(), t0, units = "secs"))) }

message("building object and scaling ...")
o0 <- CreateSeuratObject(CreateAssayObject(counts = Matrix(cnt, sparse = TRUE)), assay = "RNA")
sc <- tm(ScaleData(o0, features = rownames(o0), verbose = FALSE)); o0 <- sc$v
message(sprintf("  ScaleData %.0f s", sc$s))

build <- function(seed) {
  set.seed(seed)
  tm(FindNeighbors(o0, dims = NULL, features = foi, verbose = FALSE,
                   graph.name = c("ckknn", "cksnn")))
}
clust <- function(obj, res, seed = SEED) {
  set.seed(seed)
  tm(FindClusters(obj, resolution = res, graph.name = "cksnn",
                  random.seed = seed, verbose = FALSE))
}

## ---- graph A, then two clusterings on it at identical settings -------------------------
gA <- build(SEED); message(sprintf("  FindNeighbors (build A) %.0f s", gA$s))
c1 <- clust(gA$v, 1.0); c2 <- clust(gA$v, 1.0)
labs1 <- as.integer(c1$v$seurat_clusters); labs2 <- as.integer(c2$v$seurat_clusters)
message(sprintf("  FindClusters %.1f s per call", mean(c(c1$s, c2$s))))

## ---- graph B: identical input and seed, rebuilt from scratch ---------------------------
gB <- build(SEED); message(sprintf("  FindNeighbors (build B) %.0f s", gB$s))
c3 <- clust(gB$v, 1.0); labs3 <- as.integer(c3$v$seurat_clusters)

ari <- function(a, b) mclust::adjustedRandIndex(a, b)
same_graph  <- ari(labs1, labs2)
rebuilt     <- ari(labs1, labs3)

## ---- resolution sweep on one built graph ------------------------------------------------
sweep_s <- 0; ks <- integer(0)
for (r in RES) { z <- clust(gA$v, r); sweep_s <- sweep_s + z$s
  ks <- c(ks, length(unique(z$v$seurat_clusters))) }

naive <- length(RES) * (gA$s + mean(c(c1$s, c2$s)))
reuse <- gA$s + sweep_s

res <- data.frame(
  depth = DEPTH,
  scale_s = round(sc$s, 1), graph_s = round(mean(c(gA$s, gB$s)), 1),
  clust_s = round(mean(c(c1$s, c2$s)), 1),
  ari_same_graph = round(same_graph, 6), ari_graph_rebuilt = round(rebuilt, 6),
  n_clusters_res1 = length(unique(labs1)),
  sweep_naive_s = round(naive, 1), sweep_reuse_s = round(reuse, 1),
  speedup = round(naive / reuse, 2))
write.csv(res, file.path(OUT, sprintf("seurat-checks-d%d.csv", DEPTH)), row.names = FALSE)
write.csv(data.frame(resolution = RES, n_clusters = ks),
          file.path(OUT, sprintf("resolution-grid-d%d.csv", DEPTH)), row.names = FALSE)

message("\n==== RESULTS ====")
message(sprintf("determinism, same graph, same seed : ARI %.6f  %s", same_graph,
                if (same_graph > 0.9999) "IDENTICAL" else "NOT identical"))
message(sprintf("determinism, graph rebuilt, same seed: ARI %.6f  %s", rebuilt,
                if (rebuilt > 0.9999) "IDENTICAL" else "NOT identical"))
message(sprintf("graph %.0f s, clustering %.1f s per resolution", mean(c(gA$s, gB$s)),
                mean(c(c1$s, c2$s))))
message(sprintf("8-resolution sweep: %.0f s naive vs %.0f s reusing the graph -> %.1fx",
                naive, reuse, naive / reuse))
message("clusters by resolution: ", paste(sprintf("%.1f=%d", RES, ks), collapse = "  "))
