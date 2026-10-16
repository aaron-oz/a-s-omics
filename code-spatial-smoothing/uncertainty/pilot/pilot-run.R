## Phase 3 pilot: run every frozen comparison on the showcase fields, with naive comparators.
## Plug-in = PUBLISHED fields (common minimum of saved posterior medians) / total count, i.e. what
## a reader of the paper's maps sees. Writes pilot/results-*.csv and pilot/pilot-results.rds.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript pilot-run.R
source("compare-lib.R")
D  <- readRDS(file.path(UDIR, "design/design-regions.rds"))
T0 <- 1.25
showcase <- c("CALR-ITGA2B", "NID1-ITGB1", "ADAM12-ITGB1", "WNT4-FZD6", "RSPO1-LGR6", "DHH-PTCH1", "VEGFA-KDR")
cross <- c("CALM1-KCNQ3", "CALM1-PDE1A")
genes <- unique(gene_case(unlist(strsplit(c(showcase, cross), "-"))))
cache <- setNames(lapply(c(genes, "Cd44"), load_draws), c(genes, "Cd44"))

## plug-in fields and bin keys
tc <- readRDS(file.path(REPO, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
key <- paste(tc$x, tc$y, sep = ".")
mn <- readRDS(file.path(REPO, "code-spatial-smoothing/controls/sweep/m-est-min-matrix.rds"))
keep <- !grepl("conflicted", mn$features); Mpub <- mn$counts[keep, ]; rownames(Mpub) <- mn$features[keep]
stopifnot(all(colnames(Mpub) == key))
plug <- function(pair) if (pair == "Cd44") tc$median else Mpub[pair, ] / tc$total.count
idx <- lapply(D$regions, function(r) match(paste(r$x, r$y, sep = "."), key))
cmp <- rbind(data.table(r1 = c("box-heart-1", "box-cord-1", "box-heart-1"), r2 = c("box-heart-2", "box-cord-2", "box-liver")),
             as.data.table(t(combn(grep("^anat-", names(idx), value = TRUE), 2)))[, .(r1 = V1, r2 = V2)])

field <- function(f) if (f == "Cd44") cache[["Cd44"]] else sig_draws(f, cache)

## ---- 1. areal averages, same field
reg <- rbindlist(lapply(c(showcase, "Cd44"), function(f) {
  M <- field(f); p <- plug(f)
  rbindlist(lapply(seq_len(nrow(cmp)), function(i) {
    k1 <- idx[[cmp$r1[i]]]; k2 <- idx[[cmp$r2[i]]]
    a <- mean(p[k1]); b <- mean(p[k2]); s <- summ_contrast(area_contrast(M, k1, k2))
    data.table(field = f, r1 = cmp$r1[i], r2 = cmp$r2[i], plug_a = a, plug_b = b, plug_ratio = a / b,
               looks_diff = max(a, b) / min(a, b) >= T0,
               wilcox_p = wilcox.test(p[k1], p[k2])$p.value, as.data.table(s))
  }))
}))
fwrite(reg, "results-regional.csv")

## ---- 2. areal averages, cross field (shape-normalized), in each region
S1 <- shape_norm(field(cross[1])); S2 <- shape_norm(field(cross[2]))
p1 <- plug(cross[1]); p1 <- p1 / mean(p1); p2 <- plug(cross[2]); p2 <- p2 / mean(p2)
regx <- rbindlist(lapply(names(idx), function(nm) {
  k <- idx[[nm]]; a <- mean(p1[k]); b <- mean(p2[k])
  s <- summ_contrast(log(colMeans(S1[k, ])) - log(colMeans(S2[k, ])))
  data.table(field = paste(cross, collapse = " vs "), region = nm, plug_a = a, plug_b = b, plug_ratio = a / b,
             looks_diff = max(a, b) / min(a, b) >= T0,
             wilcox_signed_p = wilcox.test(p1[k], p2[k], paired = TRUE)$p.value, as.data.table(s))
}))
fwrite(regx, "results-crossfield-regional.csv")

## ---- 3. slices, same field: deviation from the slice average
sl <- list(); sls <- list()
for (f in c(showcase, "Cd44")) {
  M <- field(f); p <- plug(f)
  for (nm in names(D$slice_bins)) {
    sb <- D$slice_bins[[nm]]; k <- match(paste(sb$x, sb$y, sep = "."), key)
    X <- slice_dev(M, k); J <- joint_summary(X); invisible(capture.output(sc <- simconf.mc(X, alpha = 0.05)))
    d <- data.table(field = f, slice = nm, pos = sb$pos, annotation = sb$annotation,
                    plug_dev = log(p[k]) - log(mean(p[k])), med = J$med, pw_lo = J$pw_lo, pw_hi = J$pw_hi,
                    sim_lo = sc$a, sim_hi = sc$b, pw_pos = J$pw_pos, pw_neg = J$pw_neg, exc_pos = J$pos, exc_neg = J$neg)
    sl[[paste(f, nm)]] <- d
    sls[[paste(f, nm)]] <- d[, .(field = f, slice = nm, n = .N,
                                 plug_maxmin = exp(max(plug_dev) - min(plug_dev)),
                                 n_pointwise = sum(pw_pos | pw_neg),
                                 n_simband = sum(sim_lo > 0 | sim_hi < 0),
                                 n_excursion = sum(exc_pos | exc_neg))]
  }
}
fwrite(rbindlist(sl), "results-slices.csv"); fwrite(rbindlist(sls), "results-slices-summary.csv")

## ---- 4. whole fields: cross field, and Cd44 against its tissue average
t0 <- Sys.time()
Jx <- joint_summary(log(S1) - log(S2))
Xc <- log(cache[["Cd44"]]); Xc <- sweep(Xc, 2, log(colMeans(cache[["Cd44"]])), "-")
Jc <- joint_summary(Xc)
cat(sprintf("whole-field excursions: %.0f s\n", as.numeric(Sys.time() - t0, units = "secs")))
wf <- rbind(
  data.table(comparison = paste(cross, collapse = " vs "), n_bins = length(Jx$pos),
             n_pointwise_pos = sum(Jx$pw_pos), n_pointwise_neg = sum(Jx$pw_neg),
             n_excursion_pos = sum(Jx$pos), n_excursion_neg = sum(Jx$neg)),
  data.table(comparison = "Cd44 vs its tissue average", n_bins = length(Jc$pos),
             n_pointwise_pos = sum(Jc$pw_pos), n_pointwise_neg = sum(Jc$pw_neg),
             n_excursion_pos = sum(Jc$pos), n_excursion_neg = sum(Jc$neg)))
fwrite(wf, "results-wholefield.csv")
saveRDS(list(Jx = Jx, Jc = Jc, plug_x = log(p1) - log(p2), plug_c = log(tc$median) - log(mean(tc$median)),
             xy = tc[, .(x, y)]), "/var/home/aoz/data/a-s-omics-uncertainty/wholefield.rds")

## ---- console summary
cat("\nRegional, same field: looks different (plug-in r >= 1.25) vs credible (95% interval excludes 1)\n")
print(reg[, .N, by = .(looks_diff, credible, wilcox_sig = wilcox_p < 0.05)][order(-looks_diff, -credible)])
cat("\nCross field, regional:\n"); print(regx[, .(region, plug_ratio = round(plug_ratio, 2), wilcox_signed_p = signif(wilcox_signed_p, 2),
                                                 post_ratio = round(post_ratio, 2), lo = round(lo, 2), hi = round(hi, 2), credible)])
cat("\nSlices:\n"); print(rbindlist(sls))
cat("\nWhole fields:\n"); print(wf)
