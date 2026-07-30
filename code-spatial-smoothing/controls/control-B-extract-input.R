## Control B, step 1: memory-safe extraction of the bin-50 sub-region raw data.
## The full pre-fit-obj-binned-50.Rdata holds a 122M-row all.dat (~9 GB) that hung
## the box before. Here we load it once, immediately subset to the sub.bound FOV and
## the modeled genes, keep only the compact columns needed for refitting, write a
## small input object, and drop the big table. Run with a memory cap in mind.

suppressMessages({library(data.table); library(splancs)})
setwd("/var/home/aoz/Dropbox/genetics/a-s-omics")
setDTthreads(4)
out.dir <- "data-outputs/mouse-embryo-raw/controls"; dir.create(out.dir, showWarnings=FALSE, recursive=TRUE)

## genes that actually have fitted fields (== the paper's modeled set we can compare to)
GG <- readRDS(file.path(out.dir, "gene-field-matrix.rds"))
genes <- colnames(GG$G)          # title-case gene names, 912 of them
rm(GG); gc()

e <- new.env()
load("data-inputs/mouse-embryo-raw/pre-fit-obj-binned-50.Rdata", envir=e)
ad <- e$all.dat; rm(e); gc()
cat("all.dat rows:", nrow(ad), " cols:", paste(names(ad), collapse=","), "\n")

## sub.bound FOV (identical to supp-4 00-launch)
sb <- as.points(matrix(c(15000,25000,25000,15000,15000, 11500,11500,18000,18000,11500), ncol=2))
ins <- which(inout(as.matrix(ad[, .(x,y)]), sb))
sub <- ad[ins]; rm(ad); gc()
sub <- sub[feat %in% genes, .(x, y, feat, feat.count, total.count)]
setkey(sub, feat)
cat("sub-region rows:", nrow(sub), " unique feats:", length(unique(sub$feat)),
    " unique locs:", nrow(unique(sub[, .(x,y)])), "\n")

saveRDS(sub, file.path(out.dir, "bin50-subregion-raw.rds"))
cat("wrote bin50-subregion-raw.rds  (",
    round(file.info(file.path(out.dir,"bin50-subregion-raw.rds"))$size/1e6,1), "MB )\n")
