## Symmetry and width, on the log scale, of the stored marginal posteriors of expected count.
##
## Supports two numbers in the Results section "The interaction fields are plug-in estimates, and
## differ from their posterior summaries", which takes each marginal as lognormal (jensen-bias.R).
## For every gene g and bin s, from the stored summaries q0.025, q0.5, q0.975 of the expected count:
##   ratio(g, s) = (log q0.975 - log q0.5) / (log q0.5 - log q0.025)   upper over lower half-width;
##                 1 for a posterior that is symmetric on the log scale
##   sigma(g, s) = (log q0.975 - log q0.025) / (2 * 1.959964)          log-scale sd under the lognormal
## Genes: all those entering the 1,477 ligand-receptor pairs (the same set as jensen-bias.R), every
## bin of the field of view. No sampling. This replaces figures first computed on an unrecorded
## random sample of 60 features.
##
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript interval-symmetry.R
## Output: interval-symmetry-genes.csv (one row per gene), interval-symmetry-summary.csv (pooled).
suppressMessages(library(data.table))
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
pdir <- file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects")
sw   <- file.path(repo, "code-spatial-smoothing/controls/sweep")
odir <- file.path(repo, "code-spatial-smoothing/uncertainty")
mn  <- readRDS(file.path(sw, "m-est-min-matrix.rds"))
mech <- mn$features[!grepl("conflicted", mn$features)]
e <- new.env(); load(file.path(sw, "fantom-hierarchy.Robj"), envir = e); fan <- get("output", e)
fi <- match(mech, gsub("_", "-", fan$MECHANISM)); stopifnot(!anyNA(fi))
need <- unique(c(fan$LIGAND[fi], fan$RECEPTOR[fi]))
files <- list.files(pdir, "-pred-poi\\.rds$"); gname <- sub("-pred-poi\\.rds$", "", files)
gmap <- setNames(gname, toupper(gname)); stopifnot(all(need %in% names(gmap)))
cat(sprintf("pairs: %d   genes: %d\n", length(mech), length(need)))
RAT <- SIG <- NULL; rows <- vector("list", length(need))
for (k in seq_along(need)) {
  ex <- readRDS(file.path(pdir, paste0(gmap[[need[k]]], "-pred-poi.rds")))$expect
  if (is.null(RAT)) RAT <- SIG <- matrix(NA_real_, nrow(ex), length(need))
  up <- log(ex$q0.975) - log(ex$q0.5); lo <- log(ex$q0.5) - log(ex$q0.025)
  RAT[, k] <- up / lo; SIG[, k] <- (up + lo) / (2 * qnorm(0.975))
  rows[[k]] <- data.table(gene = gmap[[need[k]]], ratio_median = median(RAT[, k]), ratio_p05 = quantile(RAT[, k], .05),
                          ratio_p95 = quantile(RAT[, k], .95), sigma_median = median(SIG[, k]), n_nonfinite = sum(!is.finite(RAT[, k])))
}
G <- rbindlist(rows); fwrite(G, file.path(odir, "interval-symmetry-genes.csv"))
qs <- c(.01, .05, .25, .5, .75, .95, .99)
S <- data.table(quantile = qs, ratio_pooled = quantile(RAT, qs, na.rm = TRUE), sigma_pooled = quantile(SIG, qs, na.rm = TRUE),
                ratio_gene_median = quantile(G$ratio_median, qs), sigma_gene_median = quantile(G$sigma_median, qs))
fwrite(S, file.path(odir, "interval-symmetry-summary.csv"))
cat(sprintf("gene-bins: %d   non-finite ratios: %d\n", length(RAT), sum(!is.finite(RAT))))
cat(sprintf("share of gene-bins with ratio in [0.85, 1.12]: %.4f\n", mean(RAT >= 0.85 & RAT <= 1.12, na.rm = TRUE)))
print(S, digits = 4)
