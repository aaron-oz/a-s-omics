## Added 2026-10-29 (AOZ review of the pivot). The Wilcoxon comparator in analyze-systematic.R
## tests the PLUG-IN (smoothed, fitted) field values of the bins in two regions. That is one
## reading people make, but the more common one tests per-bin scores computed from the data
## themselves, with no spatial model. This script adds that comparator for the 13 frozen regional
## comparisons of every pair: a two-sided Wilcoxon rank-sum test on the OBSERVED per-bin
## interaction score min(y_L, y_R) / N (same common-minimum operator, applied to the observed
## normalized counts), bins treated as independent. Also records the share of bins in each region
## whose observed score is zero, since sparse partners make that score zero at most bins.
## Output: results/regional-rawwilcox.csv (pair, r1, r2, raw_wilcox_p, zero_share_1, zero_share_2,
## raw_ratio = ratio of regional means of the observed score, NA when a region mean is 0).
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript raw-wilcoxon.R
suppressMessages(library(data.table))
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"; ud <- file.path(repo, "code-spatial-smoothing/uncertainty")
I <- readRDS(file.path(ud, "klone/stage/design-index.rds"))
P <- fread(file.path(ud, "klone/stage/pairs.csv"))
tc <- readRDS(file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
key <- paste(tc$x, tc$y, sep = ".")
sub <- readRDS(file.path(repo, "data-outputs/mouse-embryo-raw/controls/bin50-subregion-raw.rds"))
sub[, k := paste(x, y, sep = ".")]
## observed rate per bin, one vector per gene, in the design bin order
genes <- unique(c(P$Lg, P$Rg))
Y <- lapply(setNames(genes, genes), function(g) { d <- sub[feat == g]; stopifnot(nrow(d) == 26000)
  r <- d$feat.count / d$total.count; r[match(key, d$k)] })
stopifnot(all(abs(tc$total.count - sub[feat == genes[1]][match(key, k), total.count]) < 1e-9))
cmp <- rbind(data.table(r1 = c("box-heart-1", "box-cord-1", "box-heart-1"), r2 = c("box-heart-2", "box-cord-2", "box-liver")),
             as.data.table(t(combn(grep("^anat-", names(I$regions), value = TRUE), 2)))[, .(r1 = V1, r2 = V2)])
out <- rbindlist(lapply(seq_len(nrow(P)), function(i) { s <- pmin(Y[[P$Lg[i]]], Y[[P$Rg[i]]])
  rbindlist(lapply(seq_len(nrow(cmp)), function(j) { k1 <- I$regions[[cmp$r1[j]]]; k2 <- I$regions[[cmp$r2[j]]]
    a <- s[k1]; b <- s[k2]
    p <- if (all(a == 0) && all(b == 0)) NA_real_ else suppressWarnings(wilcox.test(a, b)$p.value)
    data.table(pair = P$pair[i], r1 = cmp$r1[j], r2 = cmp$r2[j], raw_wilcox_p = p,
               zero_share_1 = mean(a == 0), zero_share_2 = mean(b == 0),
               raw_ratio = if (mean(b) > 0) mean(a) / mean(b) else NA_real_) })) }))
fwrite(out, file.path(ud, "klone/results/regional-rawwilcox.csv"))
cat("instances:", nrow(out), " raw Wilcoxon p < 0.05:", out[, sum(raw_wilcox_p < .05, na.rm = TRUE)],
    " NA (both regions all zero):", out[, sum(is.na(raw_wilcox_p))], "\n")
cat("median zero share, boxes:", out[grepl("box", r1), median(c(zero_share_1, zero_share_2))],
    " anatomical:", out[grepl("anat", r1), median(c(zero_share_1, zero_share_2))], "\n")
