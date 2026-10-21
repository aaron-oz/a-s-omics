## Is D = 1,000 enough? Split each field's draws into halves (1-500, 501-1000), recompute every
## regional 95% interval on each half, and compare endpoints on the log scale. The full-D Monte
## Carlo error of an endpoint is roughly the half-vs-half difference / 2 (two independent halves,
## each with sqrt(2) times the full-D error). Also counts credible-call flips between halves.
source("compare-lib.R")
D <- readRDS(file.path(UDIR, "design/design-regions.rds"))
showcase <- c("CALR-ITGA2B", "NID1-ITGB1", "ADAM12-ITGB1", "WNT4-FZD6", "RSPO1-LGR6", "DHH-PTCH1", "VEGFA-KDR")
genes <- unique(gene_case(unlist(strsplit(showcase, "-"))))
cache <- setNames(lapply(genes, load_draws), genes)
tc <- readRDS(file.path(REPO, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
key <- paste(tc$x, tc$y, sep = ".")
idx <- lapply(D$regions, function(r) match(paste(r$x, r$y, sep = "."), key))
cmp <- rbind(data.table(r1 = c("box-heart-1", "box-cord-1", "box-heart-1"), r2 = c("box-heart-2", "box-cord-2", "box-liver")),
             as.data.table(t(combn(grep("^anat-", names(idx), value = TRUE), 2)))[, .(r1 = V1, r2 = V2)])
out <- rbindlist(lapply(showcase, function(f) { M <- sig_draws(f, cache)
  rbindlist(lapply(seq_len(nrow(cmp)), function(i) {
    c <- area_contrast(M, idx[[cmp$r1[i]]], idx[[cmp$r2[i]]])
    a <- quantile(c[1:500], c(.025, .975)); b <- quantile(c[501:1000], c(.025, .975))
    data.table(field = f, d_lo = abs(a[1] - b[1]), d_hi = abs(a[2] - b[2]), width = diff(quantile(c, c(.025, .975))),
               flip = (a[1] > 0 | a[2] < 0) != (b[1] > 0 | b[2] < 0))
  })) }))
cat(sprintf("91 instances. Half-vs-half endpoint difference (log scale): median %.3f, 95th pct %.3f, max %.3f\n",
            median(c(out$d_lo, out$d_hi)), quantile(c(out$d_lo, out$d_hi), .95), max(c(out$d_lo, out$d_hi))))
cat(sprintf("Relative to interval width: median %.3f, max %.3f.  Credible-call flips between halves: %d of 91\n",
            median(c(out$d_lo, out$d_hi) / out$width), max(c(out$d_lo, out$d_hi) / out$width), sum(out$flip)))
