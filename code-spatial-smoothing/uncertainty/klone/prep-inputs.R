## Prepare small inputs for the klone systematic run (run locally).
## Writes to klone/stage/: genes.txt (genes needed by the 1,477 pairs, title case), pairs.csv
## (pair, ligand gene, receptor gene), check-bins.rds (2,000 fixed bins x saved lambda median and
## 95% quantiles per gene, for the on-cluster reproduction check), design-index.rds (bin indices of
## the frozen regions and slices in the saved bin order).
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript prep-inputs.R
suppressMessages(library(data.table))
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
dir.create("stage", showWarnings = FALSE)
pdir <- file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects")
mn <- readRDS(file.path(repo, "code-spatial-smoothing/controls/sweep/m-est-min-matrix.rds"))
mech <- mn$features[!grepl("conflicted", mn$features)]
e <- new.env(); load(file.path(repo, "code-spatial-smoothing/controls/sweep/fantom-hierarchy.Robj"), envir = e)
fan <- as.data.table(get("output", e))[, key := gsub("_", "-", MECHANISM)]
f <- fan[match(mech, key), .(pair = key, L = LIGAND, R = RECEPTOR)]
gname <- sub("-pred-poi\\.rds$", "", list.files(pdir, "-pred-poi\\.rds$")); gmap <- setNames(gname, toupper(gname))
f[, `:=`(Lg = gmap[L], Rg = gmap[R])]; stopifnot(!anyNA(f$Lg), !anyNA(f$Rg))
genes <- sort(unique(c(f$Lg, f$Rg)))
sub <- readRDS(file.path(repo, "data-outputs/mouse-embryo-raw/controls/bin50-subregion-raw.rds"))
stopifnot(all(genes %in% sub$feat))
writeLines(genes, "stage/genes.txt"); fwrite(f[, .(pair, Lg, Rg)], "stage/pairs.csv")
cat("pairs", nrow(f), " genes", length(genes), "\n")
set.seed(1); cb <- sort(sample(26000, 2000))
chk <- lapply(setNames(genes, genes), function(g) { l <- readRDS(file.path(pdir, paste0(g, "-pred-poi.rds")))$lambda
  as.matrix(l[cb, .(median, q0.025, q0.975)]) })
saveRDS(list(bins = cb, saved = chk), "stage/check-bins.rds")
D <- readRDS(file.path(repo, "code-spatial-smoothing/uncertainty/design/design-regions.rds"))
tc <- readRDS(file.path(pdir, "Cd44-pred-poi.rds"))$lambda; key <- paste(tc$x, tc$y, sep = ".")
idx <- lapply(D$regions, function(r) match(paste(r$x, r$y, sep = "."), key))
sidx <- lapply(D$slice_bins, function(s) match(paste(s$x, s$y, sep = "."), key))
saveRDS(list(regions = idx, slices = sidx), "stage/design-index.rds")
print(file.info(list.files("stage", full.names = TRUE))[, "size", drop = FALSE])
