## Added 2026-10-29 (review plan item P3): Monte Carlo stability of the whole-field and block
## excursion sets of Fig 8 (three comparisons). Each joint set is recomputed on draws 1-500 and
## 501-1000 separately and compared with the full-D set. Contrast matrices built exactly as in
## additions-run.R. Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript wholefield-stability.R
source("compare-lib.R")
tc <- readRDS(file.path(REPO, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
cache <- setNames(lapply(c("Calm1", "Kcnq3", "Pde1a", "Cd44", "Pkm"), load_draws), c("Calm1", "Kcnq3", "Pde1a", "Cd44", "Pkm"))
field <- function(f) if (f == "Cd44") cache[["Cd44"]] else sig_draws(f, cache)
blk <- (ceiling(((tc$x - 15025) / 50 + 1) / 5) - 1) * 26 + ceiling(((tc$y - 11525) / 50 + 1) / 5)
bmean <- function(X) rowsum(X, blk, reorder = TRUE) / 25
mk <- list()
S1 <- shape_norm(field("CALM1-KCNQ3")); S2 <- shape_norm(field("CALM1-PDE1A"))
mk[["CALM1-KCNQ3 vs CALM1-PDE1A"]] <- list(bin = log(S1) - log(S2), blk = log(bmean(S1)) - log(bmean(S2)))
for (f in c("Cd44", "PKM-CD44")) { S <- shape_norm(field(f)); mk[[paste(f, "vs its tissue average")]] <- list(bin = log(S), blk = log(bmean(S))) }
rm(S1, S2); gc(verbose = FALSE)
sets <- function(X) { J <- joint_summary(X); J$pos | J$neg }
out <- rbindlist(lapply(names(mk), function(nm) rbindlist(lapply(c("bin", "blk"), function(lev) {
  X <- mk[[nm]][[lev]]; full <- sets(X); h1 <- sets(X[, 1:500]); h2 <- sets(X[, 501:1000])
  data.table(comparison = nm, level = lev, n = nrow(X), full = sum(full), half1 = sum(h1), half2 = sum(h2),
             both_halves = sum(h1 & h2), half1_in_full = sum(h1 & full), half2_in_full = sum(h2 & full))
}))))
fwrite(out, "results-wholefield-stability.csv"); print(out)
