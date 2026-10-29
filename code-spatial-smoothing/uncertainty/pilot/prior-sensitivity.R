## REVIEW-2026-10-29.org item P2: the 91 showcase regional contrasts under the baseline prior and
## the three alternatives fitted by refit-draws-prior.R. Reports changed credible calls, interval
## width ratios and the largest change in posterior median ratio.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript prior-sensitivity.R
source("compare-lib.R")
D <- readRDS(file.path(UDIR, "design/design-regions.rds"))
showcase <- c("CALR-ITGA2B", "NID1-ITGB1", "ADAM12-ITGB1", "WNT4-FZD6", "RSPO1-LGR6", "DHH-PTCH1", "VEGFA-KDR")
genes <- unique(gene_case(unlist(strsplit(showcase, "-"))))
tc <- readRDS(file.path(REPO, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
key <- paste(tc$x, tc$y, sep = ".")
idx <- lapply(D$regions, function(r) match(paste(r$x, r$y, sep = "."), key))
cmp <- rbind(data.table(r1 = c("box-heart-1", "box-cord-1", "box-heart-1"), r2 = c("box-heart-2", "box-cord-2", "box-liver")),
             as.data.table(t(combn(grep("^anat-", names(idx), value = TRUE), 2)))[, .(r1 = V1, r2 = V2)])
run <- function(dir) { DDIR <<- dir; cache <- setNames(lapply(genes, load_draws), genes)
  hyp <- rbindlist(lapply(genes, function(g) { h <- readRDS(file.path(dir, paste0(g, "-lambda-draws.rds")))$hyper
    data.table(gene = g, range = h[1, "mean"], sigma = h[2, "mean"]) }))
  res <- rbindlist(lapply(showcase, function(f) { M <- sig_draws(f, cache)
    rbindlist(lapply(seq_len(nrow(cmp)), function(i) { s <- summ_contrast(area_contrast(M, idx[[cmp$r1[i]]], idx[[cmp$r2[i]]]))
      data.table(field = f, r1 = cmp$r1[i], r2 = cmp$r2[i], med = s$post_ratio, lo = s$lo, hi = s$hi, credible = s$credible) })) }))
  list(res = res, hyp = hyp) }
base <- run("/var/home/aoz/data/a-s-omics-uncertainty/draws")
labs <- c(a = "range prior halved, P(range < 250) = 0.95", b = "range prior doubled, P(range < 1000) = 0.95", c = "sigma prior relaxed, P(sigma > 1) = 0.05")
out <- rbindlist(lapply(names(labs), function(p) { alt <- run(paste0("/var/home/aoz/data/a-s-omics-uncertainty/draws-prior-", p))
  m <- merge(base$res, alt$res, by = c("field", "r1", "r2"), suffixes = c("_0", "_1"))
  h <- merge(base$hyp, alt$hyp, by = "gene", suffixes = c("_0", "_1"))
  data.table(prior = p, label = labs[p], n = nrow(m), calls_changed = sum(m$credible_0 != m$credible_1),
             to_credible = sum(!m$credible_0 & m$credible_1), to_not_credible = sum(m$credible_0 & !m$credible_1),
             width_ratio_median = median((log(m$hi_1) - log(m$lo_1)) / (log(m$hi_0) - log(m$lo_0))),
             width_ratio_min = min((log(m$hi_1) - log(m$lo_1)) / (log(m$hi_0) - log(m$lo_0))),
             width_ratio_max = max((log(m$hi_1) - log(m$lo_1)) / (log(m$hi_0) - log(m$lo_0))),
             max_abs_log_shift_median = max(abs(log(m$med_1) - log(m$med_0))),
             range_ratio_median = median(h$range_1 / h$range_0), sigma_ratio_median = median(h$sigma_1 / h$sigma_0)) }))
fwrite(out, "results-prior-sensitivity.csv"); print(out, digits = 3)
