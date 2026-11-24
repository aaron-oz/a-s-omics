## Effective number of independent bins in a region, from the joint posterior draws.
##
## For a field M with D joint posterior draws and a region of n bins, the regional average in draw d
## is A^(d) = mean over the n bins of M^(d)(s). If the bins' posteriors were independent,
## Var(A) would equal mean_s Var(M(s)) / n. We define the effective number of independent bins as
##     n_eff = mean_s Var(M(s)) / Var(A),
## so that n_eff = n under independence and n_eff = 1 if every bin moved together. Var(.) is the
## variance over draws. This describes the posterior uncertainty of the field (which is what limits
## a statement about a regional average), not the raw counts. A second version on the log scale,
## with M replaced by log M, is reported as a check, because rate-scale variances are right-skewed.
## Fields: the seven example pairs and Cd44. Regions: the five 15 x 15 boxes and five annotated tissues.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript effective-n.R
## Output: results-effective-n.csv
source("compare-lib.R")
D  <- readRDS(file.path(UDIR, "design/design-regions.rds"))
showcase <- c("CALR-ITGA2B", "NID1-ITGB1", "ADAM12-ITGB1", "WNT4-FZD6", "RSPO1-LGR6", "DHH-PTCH1", "VEGFA-KDR")
genes <- unique(gene_case(unlist(strsplit(showcase, "-"))))
cache <- setNames(lapply(c(genes, "Cd44"), load_draws), c(genes, "Cd44"))
tc <- readRDS(file.path(REPO, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
key <- paste(tc$x, tc$y, sep = ".")
idx <- lapply(D$regions, function(r) match(paste(r$x, r$y, sep = "."), key))
field <- function(f) if (f == "Cd44") cache[["Cd44"]] else sig_draws(f, cache)
neff <- function(X) mean(apply(X, 1, var)) / var(colMeans(X))
out <- rbindlist(lapply(c(showcase, "Cd44"), function(f) { M <- field(f)
  rbindlist(lapply(names(idx), function(nm) { X <- M[idx[[nm]], , drop = FALSE]
    data.table(field = f, region = nm, n_bins = nrow(X), draws = ncol(X), neff_rate = neff(X), neff_log = neff(log(X))) })) }))
fwrite(out, "results-effective-n.csv")
print(out[grepl("^box", region)], digits = 3)
bx <- out[grepl("^box", region) & field != "Cd44"]
cat(sprintf("\nBoxes (225 bins), seven example pairs, %d instances: n_eff rate scale median %.1f, range %.1f to %.1f; log scale median %.1f, range %.1f to %.1f\n",
            nrow(bx), median(bx$neff_rate), min(bx$neff_rate), max(bx$neff_rate), median(bx$neff_log), min(bx$neff_log), max(bx$neff_log)))
cd <- out[grepl("^box", region) & field == "Cd44"]
cat(sprintf("Boxes, Cd44: rate scale %.1f to %.1f; log scale %.1f to %.1f\n", min(cd$neff_rate), max(cd$neff_rate), min(cd$neff_log), max(cd$neff_log)))
an <- out[!grepl("^box", region) & field != "Cd44"]
print(an[, .(n_bins = n_bins[1], neff_rate_median = median(neff_rate), neff_log_median = median(neff_log)), by = region], digits = 3)
