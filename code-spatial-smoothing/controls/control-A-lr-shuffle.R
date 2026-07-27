## Control A (R2.8): ligand-receptor PAIR-SHUFFLE null, no model re-fit.
## Reuses the 912 fitted per-gene Poisson fields (expect$median) from the paper's
## 2025-04-10 run. Builds the common-minimum convolution for the TRUE FANTOM5
## pairings vs. permuted (scrambled) pairings, runs the paper's exact clustering
## (ScaleData -> FindNeighbors/FindClusters at resolution 1, no PCA), and scores
## spatial coherence with ELSA (raster, categorical) at d=50 and d=100.
## Output is summary-only (cluster count + mean ELSA per run): a few KB.
##
## Tests whether the TRUE L-R pairing yields more/coherent spatial domains than
## scrambled pairings that share identical field statistics, dimensionality, and
## clustering -- i.e. structure "beyond smoothing/dimensionality/clustering".
##
## usage: Rscript control-A-lr-shuffle.R <n_shuffles> <seed>

suppressMessages({
  library(data.table); library(RcppHNSW); library(igraph); library(Matrix)
  library(elsa); library(raster)
})

## Self-contained clustering matching the paper's recipe (scale features, no PCA,
## KNN graph in the full feature space -> shared-NN Jaccard -> Louvain at res 1).
## Uses RcppHNSW + igraph instead of Seurat's neighbor backend, which aborts on a
## libstdc++ assertion in this distrobox build. Applied identically to observed
## and shuffled pairings, so the comparison is valid.
cluster_labels <- function(M, k=20, res=1){
  Xs <- scale(M); Xs[is.na(Xs)] <- 0; Xs[Xs>10] <- 10; Xs[Xs< -10] <- -10  # ScaleData-style, clip at 10
  idx <- hnsw_knn(Xs, k=k, distance="euclidean", n_threads=4)$idx
  n <- nrow(Xs)
  A <- sparseMatrix(i=rep(seq_len(n),each=k), j=as.vector(t(idx)), x=1, dims=c(n,n))
  ov <- as(A %*% t(A), "TsparseMatrix")                                    # shared-neighbor counts
  jac <- ov@x/(2*k - ov@x); keep <- jac > (1/15) & ov@i != ov@j           # Jaccard, Seurat prune 1/15
  g <- simplify(graph_from_data_frame(
        data.frame(from=ov@i[keep]+1, to=ov@j[keep]+1, weight=jac[keep]), directed=FALSE))
  memb <- rep(NA_integer_, n); cl <- cluster_louvain(g, resolution=res)
  memb[as.integer(V(g)$name)] <- membership(cl)
  if(any(is.na(memb))) memb[is.na(memb)] <- max(memb,na.rm=TRUE) + seq_len(sum(is.na(memb)))
  memb
}
args <- commandArgs(trailingOnly=TRUE)
n.shuf <- if(length(args)>=1) as.integer(args[1]) else 3L
seed0  <- if(length(args)>=2) as.integer(args[2]) else 100L

setwd("/var/home/aoz/Dropbox/genetics/a-s-omics")
pred.dir <- "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects"
out.dir  <- "data-outputs/mouse-embryo-raw/controls"
dir.create(out.dir, recursive=TRUE, showWarnings=FALSE)
cache.G  <- file.path(out.dir, "gene-field-matrix.rds")

tic <- function(){ now <- proc.time(); function() round((proc.time()-now)[3],1) }

## ---- 1. assemble per-gene field matrix G (locations x genes), cached ----
if(file.exists(cache.G)){
  message("loading cached gene-field matrix"); GG <- readRDS(cache.G)
  G <- GG$G; xy <- GG$xy
} else {
  T0 <- tic()
  files <- list.files(pred.dir, pattern="-pred-poi\\.rds$", full.names=TRUE)
  genes <- sub("-pred-poi\\.rds$", "", basename(files))
  p1 <- readRDS(files[1])
  xy <- as.matrix(p1$expect[, .(x,y)])
  nloc <- nrow(xy)
  G <- matrix(NA_real_, nrow=nloc, ncol=length(files), dimnames=list(NULL, genes))
  for(i in seq_along(files)){
    p <- readRDS(files[i])
    G[, i] <- p$expect$median
    if(i %% 200 == 0) message(sprintf("  loaded %d/%d fields", i, length(files)))
  }
  saveRDS(list(G=G, xy=xy), cache.G)
  message(sprintf("built + cached G [%d x %d] in %ss", nrow(G), ncol(G), T0()))
}
genes <- colnames(G)
GENE <- toupper(genes)

## ---- 2. FANTOM5 mechanisms with both partners present in G ----
load("data-inputs/mouse-embryo-raw/fantom-hierarchy.Robj")
hier <- as.data.table(get(ls()[sapply(ls(), function(o) is.data.frame(get(o)))][1]))
mech <- unique(hier[, .(L=LIGAND, R=RECEPTOR)])
mech <- mech[L %in% GENE & R %in% GENE]
mech[, li := match(L, GENE)][, ri := match(R, GENE)]
mech <- mech[!is.na(li) & !is.na(ri)]
message(sprintf("mechanisms with both partners present: %d", nrow(mech)))

## ---- 3. convolution + clustering + ELSA evaluator ----
# rasterize once-derived helper; coords are a regular grid
elsa.mean <- function(labels, dd){
  r <- rasterFromXYZ(cbind(xy[,1], xy[,2], as.integer(labels)))
  e <- elsa(r, d=dd, categorical=TRUE)
  mean(e@data@values, na.rm=TRUE)
}

eval.pairing <- function(li, ri, tag){
  # common-minimum convolution: pmin of the two gene fields, per mechanism (loc x nmech)
  M <- pmin(G[, li, drop=FALSE], G[, ri, drop=FALSE])
  labs <- cluster_labels(M)
  data.frame(tag=tag, n_clusters=length(unique(labs)),
             elsa_d50=elsa.mean(labs,50), elsa_d100=elsa.mean(labs,100))
}

## ---- 4. observed (true pairing) + shuffled nulls ----
results <- list()
Tt <- tic()
results[["observed"]] <- eval.pairing(mech$li, mech$ri, "observed")
message(sprintf("observed: %d clusters, ELSA(d50)=%.3f  [%ss]",
                results[["observed"]]$n_clusters, results[["observed"]]$elsa_d50, Tt()))

for(b in seq_len(n.shuf)){
  set.seed(seed0 + b)
  ri.shuf <- sample(mech$ri)                 # scramble receptor assignment across mechanisms
  Ts <- tic()
  results[[paste0("shuffle", b)]] <- eval.pairing(mech$li, ri.shuf, paste0("shuffle", b))
  message(sprintf("shuffle %d/%d: %d clusters, ELSA(d50)=%.3f  [%ss]",
                  b, n.shuf, results[[paste0("shuffle",b)]]$n_clusters,
                  results[[paste0("shuffle",b)]]$elsa_d50, Ts()))
}

res <- rbindlist(results)
outfile <- file.path(out.dir, "control-A-lr-shuffle.csv")
fwrite(res, outfile)
message("wrote ", outfile)
print(res)
