## Second reproduction check, robust to Monte Carlo error that is shared across bins (which the
## stage-1 z-statistic wrongly treats as independent). Per gene, on the 2,000 fixed check bins:
## med_dlog = median over bins of (log refit median - log saved median), and width_ratio = median
## over bins of (refit 95% width / saved 95% width) on the log scale. Seed-to-seed spread measured
## locally on 7 genes (2026-10-21): |med_dlog| <= 0.019, width_ratio 0.96-1.04. Flag (soft, not an
## exclusion) when |med_dlog| > 0.04 or width_ratio outside 0.9-1.1.
suppressMessages(library(data.table))
root <- Sys.getenv("ROOT"); ddir <- Sys.getenv("DRAWS"); odir <- Sys.getenv("OUT")
chk <- readRDS(file.path(root, "check-bins.rds")); genes <- readLines(file.path(root, "genes.txt"))
res <- rbindlist(lapply(genes, function(g) { f <- file.path(ddir, paste0(g, ".i16")); if (!file.exists(f)) return(NULL)
  con <- file(f, "rb"); v <- readBin(con, "integer", 26000L * 1000L, size = 2); close(con)
  L <- matrix(v, 26000, 1000)[chk$bins, ] / 1000 - 20; sv <- chk$saved[[g]]
  q <- t(apply(L, 1, quantile, c(.025, .5, .975)))
  data.table(gene = g, med_dlog = median(q[, 2] - log(sv[, 1])),
             width_ratio = median((q[, 3] - q[, 1]) / (log(sv[, 3]) - log(sv[, 2]))))
}))
res[, flag := abs(med_dlog) > 0.04 | width_ratio < 0.9 | width_ratio > 1.1]
fwrite(res, file.path(odir, "check2.csv"))
