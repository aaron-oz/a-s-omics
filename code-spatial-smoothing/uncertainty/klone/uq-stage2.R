## Systematic run, stage 2 (one array task = one batch of ligand-receptor pairs). Reads the two
## genes' int16 draws, forms M = min(lambda_L, lambda_R) per bin and draw, and writes per pair
## (small; everything else is computed locally from these):
##   tissue mean of M per draw (1,000); mean of M over each frozen region per draw (10 x 1,000);
##   the same region means after the decorrelating shuffle (each bin's draws permuted
##   independently, seed per pair); log M at the slice bins (409 x 1,000), as int16.
## Idempotent with temp-file renames, like stage 1.
## env: ROOT, DRAWS, OUT, BATCH
suppressMessages(library(data.table))
root <- Sys.getenv("ROOT"); ddir <- Sys.getenv("DRAWS"); odir <- file.path(Sys.getenv("OUT"), "pairs")
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0")); B <- as.integer(Sys.getenv("BATCH", "25"))
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
P <- fread(file.path(root, "pairs.csv")); I <- readRDS(file.path(root, "design-index.rds"))
sl <- unlist(I$slices, use.names = FALSE)
rd <- function(g) { con <- file(file.path(ddir, paste0(g, ".i16")), "rb"); v <- readBin(con, "integer", 26000L * 1000L, size = 2); close(con)
  matrix(v / 1000 - 20, 26000, 1000) }
rows <- (task * B + 1):min(nrow(P), (task + 1) * B)
for (r in rows) {
  p <- P[r]; fo <- file.path(odir, paste0(p$pair, ".rds")); if (file.exists(fo)) next
  M <- exp(pmin(rd(p$Lg), rd(p$Rg)))
  set.seed(20261021L + r); Md <- t(apply(M, 1, sample))
  out <- list(pair = p$pair, tissue_mean = colMeans(M),
              region_mean = t(sapply(I$regions, function(k) colMeans(M[k, , drop = FALSE]))),
              region_mean_decor = t(sapply(I$regions, function(k) colMeans(Md[k, , drop = FALSE]))),
              slice_logM_i16 = matrix(as.integer(round((log(M[sl, ]) + 20) * 1000)), length(sl)))
  saveRDS(out, paste0(fo, ".tmp")); file.rename(paste0(fo, ".tmp"), fo)
}
