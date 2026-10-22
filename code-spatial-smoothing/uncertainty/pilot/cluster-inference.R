## Posterior-draw cluster inference (DESIGN.org, addition 4, logged before running): u = 0, k = 25.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript cluster-inference.R
source("compare-lib.R"); suppressMessages(library(terra))
U <- 0; K <- 25
tc <- readRDS(file.path(REPO, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
ii <- as.integer((tc$x - 15025) / 50 + 1); jj <- as.integer((tc$y - 11525) / 50 + 1)
cell <- (130L - jj) * 200L + ii                       # terra cell number (row 1 = top = largest y)
cache <- setNames(lapply(c("Calm1", "Kcnq3", "Pde1a", "Cd44", "Pkm"), load_draws), c("Calm1", "Kcnq3", "Pde1a", "Cd44", "Pkm"))
field <- function(f) if (f == "Cd44") cache[["Cd44"]] else sig_draws(f, cache)
X <- list()
X[["CALM1-KCNQ3 vs CALM1-PDE1A"]] <- log(shape_norm(field("CALM1-KCNQ3"))) - log(shape_norm(field("CALM1-PDE1A")))
for (f in c("Cd44", "PKM-CD44")) X[[paste(f, "vs its tissue average")]] <- log(shape_norm(field(f)))
r0 <- rast(nrows = 130, ncols = 200, xmin = 0, xmax = 200, ymin = 0, ymax = 130)
big_component <- function(mark) {                     # mark: logical by bin -> logical "in a component >= K bins"
  v <- rep(NA_real_, 26000); v[cell[mark]] <- 1
  p <- patches(setValues(r0, v), directions = 4)
  lab <- values(p)[cell]; sz <- table(lab)
  !is.na(lab) & as.vector(sz[as.character(lab)]) >= K
}
out <- list(); res <- list()
for (nm in names(X)) {
  x <- X[[nm]]; t0 <- Sys.time()
  Ip <- sapply(seq_len(ncol(x)), function(d) big_component(x[, d] > U))
  In <- sapply(seq_len(ncol(x)), function(d) big_component(x[, d] < -U))
  ep <- excursions.mc(Ip * 1, alpha = .05, u = .5, type = ">"); en <- excursions.mc(In * 1, alpha = .05, u = .5, type = ">")
  res[[nm]] <- list(pos = as.logical(ep$E), neg = as.logical(en$E), p_pos = rowMeans(Ip), p_neg = rowMeans(In))
  out[[nm]] <- data.table(comparison = nm, cluster_joint_pos = sum(res[[nm]]$pos), cluster_joint_neg = sum(res[[nm]]$neg),
                          secs = round(as.numeric(Sys.time() - t0, units = "secs")))
}
A <- readRDS("/var/home/aoz/data/a-s-omics-uncertainty/wholefield-additions.rds")
tab <- rbindlist(out)[, `:=`(bin_joint = sapply(comparison, function(n) sum(A$maps[[paste(n, "bin")]]$J$pos | A$maps[[paste(n, "bin")]]$J$neg)),
                            block_joint_bins = sapply(comparison, function(n) 25 * sum(A$maps[[paste(n, "blk")]]$J$pos | A$maps[[paste(n, "blk")]]$J$neg)))]
fwrite(tab, "results-cluster-inference.csv"); print(tab)
saveRDS(res, "/var/home/aoz/data/a-s-omics-uncertainty/cluster-inference.rds")
