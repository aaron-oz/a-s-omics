## Seurat replication of the random-order feature-depth sweep.
##
## The main sweep (sweep-run.R) used the RcppHNSW + igraph clustering substitute built for
## the R2.8 controls, because Seurat's FindNeighbors aborts on a libstdc++ assertion in the
## local emacs-r container. It found that every one of 99 random feature subsets was more
## spatially coherent than the published variance-ordered selection at depths 150 to 1000.
##
## That result cannot go into the manuscript until it is shown not to be an artifact of the
## substitute. This script runs the same comparison through Seurat's own neighbour and
## clustering path, mirroring `Sup7 - Clustering vs Features.R`: ScaleData over all features,
## then FindNeighbors(dims = NULL, features = foi) so the graph is built in the full feature
## space with no PCA, then FindClusters at resolution 1.
##
## Reduced design: three depths, one ordered arm, one reverse arm, and 19 random draws. The
## effect being tested is roughly a factor of two, against a run-to-run spread previously
## measured at about 0.0045 ELSA units, so 19 draws is ample to see it or rule it out.
##
##   SLURM:  SLURM_ARRAY_TASK_ID=3 Rscript sweep-seurat.R
##   local:  Rscript sweep-seurat.R 0 20

## The apptainer image puts its own site-library ahead of R_LIBS_USER, so the container's
## ggplot2 3.5.1 shadows the 4.0.3 in Rlib-ctrlB and Seurat 5.5.1 fails to load looking for
## ggplot2::is_ggplot, which only exists from 3.5.2. Prepend the user library explicitly.
local({
  ul <- Sys.getenv("R_LIBS_USER")
  if (nzchar(ul) && dir.exists(ul)) .libPaths(c(ul, .libPaths()))
})

## Seurat dispatches FindNeighbors through future, whose default 500 MiB cap on exported
## globals is exceeded at high feature depth (873 MiB at depth 1300). Sup5 raises this to
## 10 GB at the top of the file; Sup7 does not, and so only runs at all if Sup5 has already
## been sourced in the same session. Set it explicitly rather than inherit it.
options(future.globals.maxSize = 10.0 * 1e9)

suppressMessages({
  library(Seurat); library(Matrix); library(elsa); library(raster)
})
message("ggplot2 ", packageVersion("ggplot2"), ", Seurat ", packageVersion("Seurat"))

ROOT   <- Sys.getenv("ASOMICS_SWEEP_ROOT", unset = ".")
OUTDIR <- file.path(ROOT, Sys.getenv("ASOMICS_SEURAT_OUT", "results-seurat"))
## Seurat costs 699 s at depth 150 and 1702 s at depth 350, so depth 1300 would dominate the
## budget while testing the weakest part of the effect (p = 0.940 there against 1.000 at
## 150-1000 in the substitute). Depths are therefore configurable, and the random arm runs at
## 150 and 350 only; the ordered arm covers all three as a reference.
## Accept ':' as well as ',' as the separator. SLURM's --export parses its own value on
## commas, so passing "150,350" through it silently becomes "150" plus a stray token, and
## the run quietly covers one depth instead of two. Use colons when submitting.
DEPTHS <- as.numeric(strsplit(Sys.getenv("ASOMICS_SEURAT_DEPTHS", "150:350"), "[,:]")[[1]])
N_RANDOM <- 19L
SEED <- 42L
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

D   <- readRDS(file.path(ROOT, "m-est-min-matrix.rds"))
dup <- duplicated(D$counts)
cnt <- D$counts[!dup, , drop = FALSE]          # mechanisms x bins, as Seurat expects
rk  <- D$rank[!dup]
xy  <- D$xy
NF  <- nrow(cnt)
ORD <- order(rk); REV <- rev(ORD)
rownames(cnt) <- make.unique(rownames(cnt))
colnames(cnt) <- paste0("bin", seq_len(ncol(cnt)))
message(sprintf("%d distinct mechanisms x %d bins", NF, ncol(cnt)))

## Build the object once; ScaleData over all features, exactly as Sup7 does before its loop.
message("building Seurat object and scaling ...")
obj <- CreateSeuratObject(CreateAssayObject(counts = Matrix(cnt, sparse = TRUE)), assay = "RNA")
obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)

elsa_at <- function(labs, ds = c(50, 100)) {
  r <- rasterFromXYZ(data.frame(x = xy$x, y = xy$y, z = as.numeric(as.factor(labs))))
  vapply(ds, function(d) mean(elsa(r, d = d, categorical = TRUE)@data@values, na.rm = TRUE), 0)
}

evaluate <- function(cols, arm, draw, depth) {
  foi <- rownames(obj)[cols]
  t0 <- Sys.time()
  set.seed(SEED)
  o <- FindNeighbors(obj, dims = NULL, features = foi, verbose = FALSE,
                     graph.name = c("swknn", "swsnn"))
  o <- FindClusters(o, resolution = 1, graph.name = "swsnn", random.seed = SEED,
                    verbose = FALSE)
  labs <- o$seurat_clusters
  e <- elsa_at(labs)
  data.frame(arm = arm, draw = draw, depth = depth,
             domains = length(unique(labs)),
             elsa_d50 = e[1], elsa_d100 = e[2],
             seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
}

task_spec <- function(id) {
  if (id == 0L) return(list(arm = "ordered", draw = 0L))
  if (id == 1L) return(list(arm = "reverse", draw = 0L))
  list(arm = "random", draw = id - 1L)
}

run_task <- function(id) {
  s <- task_spec(id)
  out <- file.path(OUTDIR, sprintf("seurat-%s-%03d.csv", s$arm, s$draw))
  rows <- list()
  for (n in DEPTHS) {
    cols <- switch(s$arm,
      ordered = ORD[seq_len(n)],
      reverse = REV[seq_len(n)],
      random  = { set.seed(2000L + s$draw); sample.int(NF, n) })
    r <- evaluate(cols, s$arm, s$draw, n)
    message(sprintf("[task %3d] %-7s draw %2d depth %5d -> %2d domains, ELSA(d50) %.4f, %.0f s",
                    id, s$arm, s$draw, n, r$domains, r$elsa_d50, r$seconds))
    rows[[length(rows) + 1L]] <- r
    ## Write after every depth. A Seurat evaluation at high depth can run for over an hour,
    ## and the first version of this script lost two completed depths when the third failed.
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  invisible(do.call(rbind, rows))
}

sid <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
if (nzchar(sid)) {
  run_task(as.integer(sid))
} else {
  a <- commandArgs(trailingOnly = TRUE)
  first <- if (length(a) >= 1) as.integer(a[1]) else 0L
  last  <- if (length(a) >= 2) as.integer(a[2]) else (1L + N_RANDOM)
  for (i in seq(first, last)) run_task(i)
}
