## Random showcase draw, rule (a) of DESIGN.org "Showcase selection rule". PLUG-IN ONLY.
## Eligible pair for a slot: plug-in ratio r >= T = 1.25 between the slot's two areas, and the
## plug-in median rate over the union of the two areas at or above the median of that quantity
## across all 1,477 pairs (excludes near-empty, speckled fields). One pair drawn per slot,
## seed 20261015. Cross-field slot: one field pair sharing a ligand, both fields passing the
## sparsity filter on the whole field of view.
## Needs plugin-rate.rds from plugin-ratios.R (TMPOUT).
## Run: TMPOUT=<dir> OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript showcase-draw.R
suppressMessages(library(data.table))
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
D <- readRDS("design-regions.rds")
P <- readRDS(file.path(Sys.getenv("TMPOUT", "."), "plugin-rate.rds"))
R <- P$R; key <- P$key; T0 <- 1.25
bins <- function(r) match(paste(r$x, r$y, sep = "."), key)
sb <- D$slice_bins[["slice-heart-to-liver"]]
areas <- list(
  "anat-heart vs anat-liver" = list(bins(D$regions[["anat-heart"]]), bins(D$regions[["anat-liver"]])),
  "box-cord-1 vs box-cord-2" = list(bins(D$regions[["box-cord-1"]]), bins(D$regions[["box-cord-2"]])),
  "slice-heart-to-liver: heart part vs liver part" = list(bins(sb[annotation == "Heart"]), bins(sb[annotation == "Liver"])))
set.seed(20261015)
picks <- rbindlist(lapply(names(areas), function(nm) {
  k1 <- areas[[nm]][[1]]; k2 <- areas[[nm]][[2]]
  a <- rowMeans(R[, k1]); b <- rowMeans(R[, k2]); r <- pmax(a, b) / pmin(a, b)
  med <- apply(R[, c(k1, k2)], 1, median)
  ok <- which(r >= T0 & med >= median(med))
  p <- ok[sample(length(ok), 1)]
  data.table(slot = nm, pair = rownames(R)[p], r = round(r[p], 2), n_eligible = length(ok))
}))
e <- new.env(); load(file.path(repo, "code-spatial-smoothing/controls/sweep/fantom-hierarchy.Robj"), envir = e)
fan <- as.data.table(get("output", e))[, key := gsub("_", "-", MECHANISM)]
lig <- fan[match(rownames(R), key), LIGAND]
medall <- apply(R, 1, median); dense <- medall >= median(medall)
cand <- rbindlist(lapply(split(which(dense), lig[dense]), function(k)
  if (length(k) >= 2) as.data.table(t(combn(k, 2))) else NULL))
cp <- cand[sample(.N, 1)]
picks <- rbind(picks, data.table(slot = "cross-field, shared ligand",
                                 pair = paste(rownames(R)[cp$V1], rownames(R)[cp$V2], sep = " & "),
                                 r = NA_real_, n_eligible = nrow(cand)))
print(picks)
fwrite(picks, "showcase-random.csv")
