## Control A replication through Seurat's own clustering path.
##
## control-A-v2.R found a clean null using the RcppHNSW + igraph substitute: the true
## ligand-receptor pairing produced no more spatially coherent domains than 99 exclusion-aware
## scrambles. The sweep replication then showed the substitute exaggerates differences between
## feature configurations by roughly tenfold relative to Seurat, which makes a null in the
## substitute conservative rather than suspect. This measures it instead of arguing it.
##
## Task 0 is the true pairing, and doubles as a second useful number: it is a full-depth
## Seurat evaluation of the paper's exact configuration, directly comparable to Control B's
## published observed value of 20 domains at ELSA(d50) 0.075, which was also produced by the
## substitute.
##
##   SLURM:  SLURM_ARRAY_TASK_ID=0 Rscript control-A-seurat.R
##   local:  Rscript control-A-seurat.R 0 3

local({
  ul <- Sys.getenv("R_LIBS_USER")
  if (nzchar(ul) && dir.exists(ul)) .libPaths(c(ul, .libPaths()))
})
options(future.globals.maxSize = 10.0 * 1e9)

suppressMessages({
  library(data.table); library(Seurat); library(Matrix); library(elsa); library(raster)
})

ROOT   <- Sys.getenv("ASOMICS_CTRLA_ROOT", unset = ".")
OUTDIR <- file.path(ROOT, Sys.getenv("ASOMICS_CTRLA_OUT", "results-ctrlA-seurat"))
SEED   <- 42L
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

GG <- readRDS(file.path(ROOT, "gene-field-matrix.rds"))
G  <- GG$G; xy <- GG$xy
GENE <- toupper(colnames(G))

e <- new.env(); load(file.path(ROOT, "fantom-hierarchy.Robj"), envir = e)
hier <- as.data.table(get(ls(e)[1], envir = e))
mech <- unique(hier[, .(L = LIGAND, R = RECEPTOR)])
mech <- mech[L %in% GENE & R %in% GENE]
mech[, li := match(L, GENE)][, ri := match(R, GENE)]
mech <- mech[!is.na(li) & !is.na(ri)]
REAL <- unique(paste(hier$LIGAND, hier$RECEPTOR, sep = "|"))
message(sprintf("%d mechanisms fitted; %d real FANTOM5 pairs", nrow(mech), length(REAL)))

## identical scramble to control-A-v2.R, so the two runs are comparable draw for draw
scramble <- function(li, ri, seed) {
  set.seed(seed)
  n <- length(li)
  ri <- sample(ri)
  is_bad <- function(a, b) (paste(GENE[a], GENE[b], sep = "|") %in% REAL) | (a == b)
  for (sw in 1:200) {
    bad <- which(is_bad(li, ri))
    if (!length(bad)) break
    for (i in bad) for (attempt in 1:100) {
      j <- sample.int(n, 1L)
      if (j == i) next
      if (!is_bad(li[i], ri[j]) && !is_bad(li[j], ri[i])) {
        tmp <- ri[i]; ri[i] <- ri[j]; ri[j] <- tmp; break
      }
    }
  }
  list(ri = ri, remaining_real = sum(is_bad(li, ri)))
}

elsa_at <- function(labs, ds = c(50, 100)) {
  r <- rasterFromXYZ(data.frame(x = xy[, 1], y = xy[, 2], z = as.numeric(as.factor(labs))))
  vapply(ds, function(d) mean(elsa(r, d = d, categorical = TRUE)@data@values, na.rm = TRUE), 0)
}

evaluate <- function(li, ri, tag, remaining = 0L) {
  t0 <- Sys.time()
  M <- pmin(G[, li, drop = FALSE], G[, ri, drop = FALSE])   # locations x mechanisms
  ndist <- sum(!duplicated(t(M)))
  cnt <- t(M)                                                # mechanisms x locations for Seurat
  rownames(cnt) <- make.unique(paste(GENE[li], GENE[ri], sep = "-"))
  colnames(cnt) <- paste0("bin", seq_len(ncol(cnt)))
  o <- CreateSeuratObject(CreateAssayObject(counts = Matrix(cnt, sparse = TRUE)), assay = "RNA")
  o <- ScaleData(o, features = rownames(o), verbose = FALSE)
  set.seed(SEED)
  o <- FindNeighbors(o, dims = NULL, features = rownames(o), verbose = FALSE,
                     graph.name = c("caknn", "casnn"))
  o <- FindClusters(o, resolution = 1, graph.name = "casnn", random.seed = SEED, verbose = FALSE)
  labs <- o$seurat_clusters
  el <- elsa_at(labs)
  data.frame(tag = tag, n_mech = length(li), n_distinct = ndist,
             remaining_real_pairs = remaining,
             n_clusters = length(unique(labs)),
             elsa_d50 = el[1], elsa_d100 = el[2],
             seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
}

run_task <- function(id) {
  if (id == 0L) {
    res <- evaluate(mech$li, mech$ri, "observed", 0L)
    out <- file.path(OUTDIR, "ctrlA-seurat-observed.csv")
  } else {
    s <- scramble(mech$li, mech$ri, seed = 5000L + id)   # same seeds as control-A-v2.R
    res <- evaluate(mech$li, s$ri, sprintf("scramble%03d", id), s$remaining_real)
    out <- file.path(OUTDIR, sprintf("ctrlA-seurat-scramble-%03d.csv", id))
  }
  write.csv(res, out, row.names = FALSE)
  message(sprintf("[%s] %d distinct, %d real left, %d clusters, ELSA(d50) %.4f, %.0f s",
                  res$tag, res$n_distinct, res$remaining_real_pairs, res$n_clusters,
                  res$elsa_d50, res$seconds))
  invisible(res)
}

sid <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
if (nzchar(sid)) run_task(as.integer(sid)) else {
  a <- commandArgs(trailingOnly = TRUE)
  for (i in seq(if (length(a) >= 1) as.integer(a[1]) else 0L,
                if (length(a) >= 2) as.integer(a[2]) else 19L)) run_task(i)
}
