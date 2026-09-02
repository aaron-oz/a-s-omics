## One-time preparation for the held-out-gene instrument.
##
## Builds a genes x bins count matrix over the paper's 26,000-bin field of view, restricted to
## genes the pipeline never used. Reads only the needed columns from the 3.8 GB whole-section
## DGE rather than loading it all, since dense it would be about 17 GB.
##
## Held-out set: the pipeline's own filter from Sup1 (expressed in at least 0.25% of bins, at
## least one bin with two counts) gives 17,120 genes; remove the 912 fitted ligand and
## receptor genes from that.

suppressMessages({library(data.table); library(Matrix)})
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
out  <- file.path(repo, "code-spatial-smoothing/controls/sweep")
dge  <- file.path(repo, "Supplements for Publication/Sup1 - Raw-to-Bin/bin.50.dge.csv")

## FOV bins, from the same matrix the clustering uses, so bin identity matches exactly
D  <- readRDS(file.path(out, "m-est-min-matrix.rds"))
xy <- D$xy
want <- paste(xy$x, xy$y, sep = ".")
message(sprintf("field of view: %d bins", length(want)))

## The column names are coordinate strings like "10025.11425", which parse as numbers, so
## fread's header auto-detection treats the header row as data. Sup1 passes header = TRUE for
## the same reason when it reads this file back.
hdr <- names(fread(dge, nrows = 0, header = TRUE))
message(sprintf("DGE columns: %d", length(hdr)))
keep <- intersect(hdr, want)
message(sprintf("matched %d of %d FOV bins in the DGE header", length(keep), length(want)))
if (length(keep) < 0.9 * length(want))
  stop("column-name convention does not match; inspect head(hdr) before proceeding")

t0 <- Sys.time()
message("reading the matched columns ...")
M <- fread(dge, select = keep, header = TRUE, showProgress = FALSE)
message(sprintf("  read %d x %d in %.0f s", nrow(M), ncol(M),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))

## Sup1 wrote the DGE without row names, so gene identity comes from the same source it used
gene_file <- file.path(repo, "data-inputs/mouse-embryo-raw/msbr-2025-04-24/split.by.gene.Robj")
e <- new.env(); load(gene_file, envir = e); genes <- names(get(ls(e)[1], envir = e))
message(sprintf("gene names from split.by.gene.Robj: %d (matrix has %d rows)",
                length(genes), nrow(M)))
stopifnot(length(genes) == nrow(M))

lr <- toupper(colnames(readRDS(file.path(repo,
        "data-outputs/mouse-embryo-raw/controls/gene-field-matrix.rds"))$G))
heldout <- which(!(toupper(genes) %in% lr))
message(sprintf("held-out genes (excluding the %d fitted): %d", length(lr), length(heldout)))

Msub <- as.matrix(M[heldout, ])
rownames(Msub) <- genes[heldout]
Msub[is.na(Msub)] <- 0
sp <- Matrix(Msub, sparse = TRUE)
saveRDS(list(counts = sp, bins = keep, genes = genes[heldout]),
        file.path(out, "heldout-gene-matrix.rds"))
message(sprintf("wrote heldout-gene-matrix.rds: %d genes x %d bins, %.1f%% nonzero, %.1f min total",
                nrow(sp), ncol(sp), 100 * Matrix::nnzero(sp) / length(sp),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
