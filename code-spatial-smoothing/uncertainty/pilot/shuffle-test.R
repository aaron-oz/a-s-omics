## Decorrelated-posterior comparator (DESIGN.org, Naive comparators, item 0).
## For each field, permute each bin's D draws independently (seed 20261016): every bin keeps its
## marginal posterior exactly; the correlation between bins is destroyed. Recompute every regional
## contrast and compare with the joint (correlated) posterior interval.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript shuffle-test.R
source("compare-lib.R")
D  <- readRDS(file.path(UDIR, "design/design-regions.rds"))
showcase <- c("CALR-ITGA2B", "NID1-ITGB1", "ADAM12-ITGB1", "WNT4-FZD6", "RSPO1-LGR6", "DHH-PTCH1", "VEGFA-KDR")
cross <- c("CALM1-KCNQ3", "CALM1-PDE1A")
genes <- unique(gene_case(unlist(strsplit(c(showcase, cross), "-"))))
cache <- setNames(lapply(c(genes, "Cd44"), load_draws), c(genes, "Cd44"))
tc <- readRDS(file.path(REPO, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
key <- paste(tc$x, tc$y, sep = ".")
idx <- lapply(D$regions, function(r) match(paste(r$x, r$y, sep = "."), key))
cmp <- rbind(data.table(r1 = c("box-heart-1", "box-cord-1", "box-heart-1"), r2 = c("box-heart-2", "box-cord-2", "box-liver")),
             as.data.table(t(combn(grep("^anat-", names(idx), value = TRUE), 2)))[, .(r1 = V1, r2 = V2)])
field <- function(f) if (f == "Cd44") cache[["Cd44"]] else sig_draws(f, cache)
set.seed(20261016)
decor <- function(M) t(apply(M, 1, sample))   # independent permutation of each bin's draws
iw <- function(s) log(s$hi / s$lo)              # interval width on the log scale

reg <- fread("results-regional.csv")
out <- rbindlist(lapply(c(showcase, "Cd44"), function(f) {
  Ms <- decor(field(f))
  rbindlist(lapply(seq_len(nrow(cmp)), function(i) {
    s <- summ_contrast(area_contrast(Ms, idx[[cmp$r1[i]]], idx[[cmp$r2[i]]]))
    data.table(field = f, r1 = cmp$r1[i], r2 = cmp$r2[i], dec_lo = s$lo, dec_hi = s$hi, dec_credible = s$credible)
  }))
}))
out <- merge(reg, out, by = c("field", "r1", "r2"))
out[, width_ratio := log(hi / lo) / log(dec_hi / dec_lo)]   # joint width / decorrelated width
fwrite(out, "results-shuffle-regional.csv")

## cross field
S1 <- shape_norm(field(cross[1])); S2 <- shape_norm(field(cross[2]))
S1d <- decor(S1); S2d <- decor(S2)
regx <- fread("results-crossfield-regional.csv")
ox <- rbindlist(lapply(names(idx), function(nm) {
  k <- idx[[nm]]; s <- summ_contrast(log(colMeans(S1d[k, ])) - log(colMeans(S2d[k, ])))
  data.table(region = nm, dec_lo = s$lo, dec_hi = s$hi, dec_credible = s$credible)
}))
ox <- merge(regx, ox, by = "region")[, width_ratio := log(hi / lo) / log(dec_hi / dec_lo)]
fwrite(ox, "results-shuffle-crossfield.csv")

sig <- out[field != "Cd44"]
cat("Signaling fields, 91 regional instances:\n")
cat(sprintf("  credible with joint posterior: %d   with decorrelated posterior: %d   Wilcoxon p<0.05: %d\n",
            sum(sig$credible), sum(sig$dec_credible), sum(sig$wilcox_p < 0.05)))
print(sig[, .N, by = .(joint = credible, decorrelated = dec_credible)])
cat("  joint / decorrelated interval width (log scale): median", round(median(sig$width_ratio), 2),
    " range", paste(round(range(sig$width_ratio), 2), collapse = " to "), "\n")
cat("  by region type: boxes median", round(median(sig[grepl("box", r1), width_ratio]), 2),
    "; anatomical median", round(median(sig[grepl("anat", r1), width_ratio]), 2), "\n")
cat("Cd44: credible joint", sum(out[field == "Cd44", credible]), " decorrelated", sum(out[field == "Cd44", dec_credible]),
    " width ratio median", round(median(out[field == "Cd44", width_ratio]), 2), "\n")
cat("Cross field, 10 regions: credible joint", sum(ox$credible), " decorrelated", sum(ox$dec_credible),
    " width ratio median", round(median(ox$width_ratio), 2), "\n")
print(sig[credible == FALSE & dec_credible == TRUE, .(field, r1, r2, plug = round(plug_ratio, 2),
      joint = sprintf("%.2f-%.2f", lo, hi), decor = sprintf("%.2f-%.2f", dec_lo, dec_hi))])
