## DESIGN.org additions logged 2026-10-21 (run after logging): coarse-grid (5 x 5 block) excursion
## sets, PKM-CD44 (drawn CD44-containing pair) vs its tissue average, and per-position
## "looks different" labels cross-tabulated against the joint excursion sets.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript additions-run.R
source("compare-lib.R")
TS <- c(1.25, 1.5, 2)
tc <- readRDS(file.path(REPO, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
key <- paste(tc$x, tc$y, sep = ".")
mn <- readRDS(file.path(REPO, "code-spatial-smoothing/controls/sweep/m-est-min-matrix.rds"))
keep <- !grepl("conflicted", mn$features); Mpub <- mn$counts[keep, ]; rownames(Mpub) <- mn$features[keep]
plug <- function(f) if (f == "Cd44") tc$median else Mpub[f, ] / tc$total.count
cache <- setNames(lapply(c("Calm1", "Kcnq3", "Pde1a", "Cd44", "Pkm"), load_draws), c("Calm1", "Kcnq3", "Pde1a", "Cd44", "Pkm"))
field <- function(f) if (f == "Cd44") cache[["Cd44"]] else sig_draws(f, cache)
blk <- (ceiling(((tc$x - 15025) / 50 + 1) / 5) - 1) * 26 + ceiling(((tc$y - 11525) / 50 + 1) / 5)   # 1..1040
bmean <- function(X) rowsum(X, blk, reorder = TRUE) / 25

## contrast draws and plug-in contrasts, bins and blocks
mk <- list()
S1 <- shape_norm(field("CALM1-KCNQ3")); S2 <- shape_norm(field("CALM1-PDE1A"))
p1 <- plug("CALM1-KCNQ3"); p1 <- p1 / mean(p1); p2 <- plug("CALM1-PDE1A"); p2 <- p2 / mean(p2)
mk[["CALM1-KCNQ3 vs CALM1-PDE1A"]] <- list(bin = log(S1) - log(S2), blk = log(bmean(S1)) - log(bmean(S2)),
  pbin = log(p1) - log(p2), pblk = log(rowsum(p1, blk) / 25) - log(rowsum(p2, blk) / 25))
for (f in c("Cd44", "PKM-CD44")) {
  M <- field(f); S <- shape_norm(M); p <- plug(f); p <- p / mean(p)
  mk[[paste(f, "vs its tissue average")]] <- list(bin = log(S), blk = log(bmean(S)), pbin = log(p), pblk = log(rowsum(p, blk) / 25))
}
rm(S1, S2); gc(verbose = FALSE)

tab <- list(); maps <- list()
for (nm in names(mk)) for (lev in c("bin", "blk")) {
  X <- mk[[nm]][[lev]]; pl <- as.vector(mk[[nm]][[paste0("p", lev)]]); J <- joint_summary(X)
  inset <- J$pos | J$neg
  maps[[paste(nm, lev)]] <- list(J = J, plug = pl)
  for (T0 in TS) {
    ld <- abs(pl) >= log(T0)
    tab[[paste(nm, lev, T0)]] <- data.table(comparison = nm, level = ifelse(lev == "bin", "bin", "5x5 block"), T = T0, n = length(pl),
      looks_diff = sum(ld), looks_diff_in_joint = sum(ld & inset), looks_diff_not_joint = sum(ld & !inset),
      looks_alike_in_joint = sum(!ld & inset), pointwise = sum(J$pw_pos | J$pw_neg), joint = sum(inset),
      looks_alike_in_joint_max_post_ratio = if (any(!ld & inset)) round(exp(max(abs(J$med[!ld & inset]))), 2) else NA_real_)
  }
}
wf <- rbindlist(tab); fwrite(wf, "results-wholefield-additions.csv")
saveRDS(list(maps = maps, xy = tc[, .(x, y)], blk = blk), "/var/home/aoz/data/a-s-omics-uncertainty/wholefield-additions.rds")
print(wf[T == 1.25])

## slices: per-position labels from the existing results (plug_dev vs the joint excursion sets)
s <- fread("results-slices.csv")
st <- rbindlist(lapply(TS, function(T0) s[, .(T = T0, n = .N, looks_diff = sum(abs(plug_dev) >= log(T0)),
  looks_diff_in_joint = sum(abs(plug_dev) >= log(T0) & (exc_pos | exc_neg)),
  looks_alike_in_joint = sum(abs(plug_dev) < log(T0) & (exc_pos | exc_neg))), by = .(field)]))
fwrite(st, "results-slices-lookslike.csv")
cat("\nSlices, all four slices pooled per field, T = 1.25:\n"); print(st[T == 1.25])
