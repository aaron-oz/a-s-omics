## AOZ option (a) for the 17 pairs touching a gene whose refit differs from the published fit:
## recompute "looks different" from the REFIT's own per-bin medians (plug-in = min of the two
## genes' refit medians), so the plug-in reading and the posterior come from the same fit.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript flagged-refit-plugin.R
suppressMessages(library(data.table))
dd <- "/var/home/aoz/data/a-s-omics-uncertainty/flagged-draws"
I <- readRDS("stage/design-index.rds"); reg <- fread("results/regional.csv")[flagged == TRUE]
P <- fread("stage/pairs.csv")[pair %in% reg$pair]
med <- function(g) { con <- file(file.path(dd, paste0(g, ".i16")), "rb"); v <- readBin(con, "integer", 26000L * 1000L, size = 2); close(con)
  exp(apply(matrix(v, 26000, 1000) / 1000 - 20, 1, median)) }
gm <- lapply(setNames(nm = unique(c(P$Lg, P$Rg))), med)
out <- rbindlist(lapply(seq_len(nrow(P)), function(i) { p <- P[i]; m <- pmin(gm[[p$Lg]], gm[[p$Rg]])
  reg[pair == p$pair][, refit_plug_ratio := mapply(function(a, b) mean(m[I$regions[[a]]]) / mean(m[I$regions[[b]]]), r1, r2)] }))
f <- function(r, T0) pmax(r, 1 / r) >= T0
out[, `:=`(ld_pub = f(plug_ratio, 1.25), ld_refit = f(refit_plug_ratio, 1.25))]
fwrite(out[, .(pair, r1, r2, plug_ratio, refit_plug_ratio, ld_pub, ld_refit, credible, lo, hi)], "results/flagged-refit-plugin.csv")
cat("flagged instances:", nrow(out), "\n")
cat("looks different (T=1.25): published fit", sum(out$ld_pub), " refit", sum(out$ld_refit), " agree", sum(out$ld_pub == out$ld_refit), "\n")
cat("looks different but not credible: published", sum(out$ld_pub & !out$credible), " refit", sum(out$ld_refit & !out$credible), "\n")
cat("looks alike but credible: published", sum(!out$ld_pub & out$credible), " refit", sum(!out$ld_refit & out$credible), "\n")
cat("median |log(refit ratio / published ratio)|:", round(median(abs(log(out$refit_plug_ratio / out$plug_ratio))), 3), "\n")
