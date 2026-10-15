## Plug-in regional ratios for the "looks different" threshold (Phase 2). PLUG-IN ONLY:
## uses the published signaling fields (min of posterior medians), converted to the rate scale
## by dividing by total count. No posterior spread is read here.
## Ratio r = max(a, b) / min(a, b), with a, b the averages of the plug-in rate field over the
## two regions' bins.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript plugin-ratios.R
suppressMessages(library(data.table))
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
D <- readRDS("design-regions.rds")
mn <- readRDS(file.path(repo, "code-spatial-smoothing/controls/sweep/m-est-min-matrix.rds"))
keep <- !grepl("conflicted", mn$features)
M <- mn$counts[keep, ]; rownames(M) <- mn$features[keep]
tc <- readRDS(file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
key <- paste(tc$x, tc$y, sep = ".")
stopifnot(all(colnames(M) == key))   # same bin order
R <- sweep(M, 2, tc$total.count, "/")   # rate-scale signaling fields, pairs x bins
cat("pairs:", nrow(R), " bins:", ncol(R), "\n")
idx <- lapply(D$regions, function(r) match(paste(r$x, r$y, sep = "."), key))
avg <- sapply(idx, function(k) rowMeans(R[, k, drop = FALSE]))   # pairs x regions
cmp <- rbind(data.table(r1 = c("box-heart-1", "box-cord-1", "box-heart-1"),
                        r2 = c("box-heart-2", "box-cord-2", "box-liver")),
             as.data.table(t(combn(grep("^anat-", names(idx), value = TRUE), 2)))[, .(r1 = V1, r2 = V2)])
out <- rbindlist(lapply(seq_len(nrow(cmp)), function(k) {
  a <- avg[, cmp$r1[k]]; b <- avg[, cmp$r2[k]]
  data.table(pair = rownames(R), r1 = cmp$r1[k], r2 = cmp$r2[k], a = a, b = b, ratio = pmax(a, b) / pmin(a, b))
}))
fwrite(out, "plugin-ratios.csv")
cat("Quantiles of the plug-in ratio, by region comparison:\n")
print(out[, as.list(round(quantile(ratio, c(.1, .25, .5, .75, .9)), 2)), by = .(r1, r2)])
cat("Share of (pair, comparison) instances with ratio >= 1.25 / 1.5 / 2:\n")
print(out[, .(ge1.25 = mean(ratio >= 1.25), ge1.5 = mean(ratio >= 1.5), ge2 = mean(ratio >= 2)), by = .(r1, r2)])
saveRDS(list(R = R, key = key, idx = idx), file.path(Sys.getenv("TMPOUT", "."), "plugin-rate.rds"))
