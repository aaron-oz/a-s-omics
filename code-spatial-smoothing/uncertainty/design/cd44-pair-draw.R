## Draw of one CD44-containing signaling pair (DESIGN.org, additions logged 2026-10-21, item 2).
## Plug-in only. Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript cd44-pair-draw.R
suppressMessages(library(data.table))
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
mn <- readRDS(file.path(repo, "code-spatial-smoothing/controls/sweep/m-est-min-matrix.rds"))
keep <- !grepl("conflicted", mn$features); M <- mn$counts[keep, ]; rownames(M) <- mn$features[keep]
tc <- readRDS(file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
R <- sweep(M, 2, tc$total.count, "/")
med <- apply(R, 1, median)
cd <- grep("(^|-)CD44(-|$)", rownames(R), value = TRUE)
ok <- cd[med[cd] >= median(med)]
cat("CD44 pairs:", length(cd), "  passing the sparsity filter:", length(ok), "\n"); print(ok)
set.seed(20261021)
pool <- if (length(ok)) ok else cd
pick <- pool[sample(length(pool), 1)]
cat("drawn:", pick, "\n")
fwrite(data.table(pair = pick, pool_size = length(pool), n_cd44_pairs = length(cd)), "cd44-pair-draw.csv")
