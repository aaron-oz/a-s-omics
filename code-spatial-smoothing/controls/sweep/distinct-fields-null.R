## Larger null for the distinct-field count of Control C (ligand-receptor pair shuffle).
##
## control-A-v2.R ran 99 shuffles, each with a full clustering. The distinct-field count needs
## no clustering: it is the number of numerically distinct columns of min(L, R) over the 1,477
## pairs. So this script reruns only that statistic, for many more shuffles, using the
## scramble() function and seeds of control-A-v2.R unchanged (seed = 5000 + id), so shuffles
## 1-99 must reproduce results-ctrlA exactly. That is checked before anything is reported.
##
## Distinctness is tested by an exact fingerprint per column (its sum and its dot product with
## a fixed random vector, compared exactly): identical columns give identical fingerprints by
## construction, and the count is validated against duplicated() on the true pairing and on
## the first 99 shuffles.
##
## Run: Rscript distinct-fields-null.R <n_shuffles> <n_cores>
suppressMessages({library(data.table); library(parallel)})
a <- commandArgs(trailingOnly = TRUE)
NS <- if (length(a) >= 1) as.integer(a[1]) else 9999L
NC <- if (length(a) >= 2) as.integer(a[2]) else 4L

GG <- readRDS("gene-field-matrix.rds"); G <- GG$G; GENE <- toupper(colnames(G))
e <- new.env(); load("fantom-hierarchy.Robj", envir = e)
hier <- as.data.table(get(ls(e)[1], envir = e))
mech <- unique(hier[, .(L = LIGAND, R = RECEPTOR)])
mech <- mech[L %in% GENE & R %in% GENE]
mech[, li := match(L, GENE)][, ri := match(R, GENE)]
mech <- mech[!is.na(li) & !is.na(ri)]
REAL <- unique(paste(hier$LIGAND, hier$RECEPTOR, sep = "|"))

## copied verbatim from control-A-v2.R
scramble <- function(li, ri, seed) {
  set.seed(seed)
  n <- length(li)
  ri <- sample(ri)
  is_bad <- function(a, b) (paste(GENE[a], GENE[b], sep = "|") %in% REAL) | (a == b)
  for (sweep in 1:200) {
    bad <- which(is_bad(li, ri))
    if (!length(bad)) break
    for (i in bad) {
      for (attempt in 1:100) {
        j <- sample.int(n, 1L)
        if (j == i) next
        if (!is_bad(li[i], ri[j]) && !is_bad(li[j], ri[i])) {
          tmp <- ri[i]; ri[i] <- ri[j]; ri[j] <- tmp; break
        }
      }
    }
  }
  list(ri = ri, remaining_real = sum(is_bad(li, ri)))
}

set.seed(123); PROJ <- runif(nrow(G))
n_distinct_fp <- function(li, ri) {
  M <- pmin(G[, li, drop = FALSE], G[, ri, drop = FALSE])
  fp <- cbind(colSums(M), as.vector(crossprod(M, PROJ)))
  sum(!duplicated(fp))
}
n_distinct_exact <- function(li, ri) sum(!duplicated(t(pmin(G[, li, drop = FALSE], G[, ri, drop = FALSE]))))

obs <- n_distinct_fp(mech$li, mech$ri)
stopifnot(obs == n_distinct_exact(mech$li, mech$ri))
cat(sprintf("pairs %d; observed distinct fields %d\n", nrow(mech), obs))

t0 <- Sys.time()
res <- rbindlist(mclapply(seq_len(NS), function(id) {
  s <- scramble(mech$li, mech$ri, seed = 5000L + id)
  data.table(id = id, n_distinct = n_distinct_fp(mech$li, s$ri), remaining_real = s$remaining_real)
}, mc.cores = NC))
cat(sprintf("%d shuffles in %.0f s on %d cores\n", NS, as.numeric(Sys.time() - t0, units = "secs"), NC))

## check against the 99 clustered shuffles of control-A-v2.R
old <- rbindlist(lapply(sprintf("results-ctrlA/ctrlA-scramble-%03d.csv", 1:min(99, NS)), fread))
chk <- merge(res, data.table(id = 1:nrow(old), old_n = old$n_distinct), by = "id")
cat(sprintf("reproduces control-A-v2 on %d of %d shuffles\n", sum(chk$n_distinct == chk$old_n), nrow(chk)))
stopifnot(all(chk$n_distinct == chk$old_n))

fwrite(res, sprintf("distinct-fields-null-%d.csv", NS))
cat(sprintf("remaining real pairs, max over shuffles: %d\n", max(res$remaining_real)))
cat(sprintf("null: mean %.2f, sd %.2f, range %d to %d\n", mean(res$n_distinct), sd(res$n_distinct),
            min(res$n_distinct), max(res$n_distinct)))
k <- sum(res$n_distinct >= obs)
cat(sprintf("shuffles >= observed: %d of %d; permutation p = (k+1)/(n+1) = %.5f; z = %.2f\n",
            k, NS, (k + 1) / (NS + 1), (obs - mean(res$n_distinct)) / sd(res$n_distinct)))
