## The calibration compares a surgery to the baseline on the SAME folds, so the quantity that
## matters is the variability of the paired DELTA across fold assignments, not of the absolute
## score. Fold noise largely cancels in a paired difference.
S <- "/tmp/claude-1000/-var-home-aoz-Dropbox-genetics-a-s-omics/49644765-1ca5-4f64-867d-f5bd061fc895/scratchpad/cvcore.R"
source(S)
ROOT <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
suppressMessages(library(data.table))
e <- new.env()
load(file.path(ROOT, "sam_sandbox/clustering vs features/m.est.min",
               "clustering vs features experiment",
               "meta.data.feature.depth.experiment.2025-07-12.Robj"), envir = e)
md <- get(ls(e)[1], envir = e)
cc <- grep("featureClusters$", colnames(md), value = TRUE)
base <- as.integer(as.factor(md[[cc[which(as.integer(sub("featureClusters$","",cc)) == 1481)]]]))

ann <- as.data.frame(fread(file.path(ROOT,
        "data-outputs/mouse-embryo-raw/controls/e16-fov-annotation.csv")))
ord <- match(paste(xy$x, xy$y, sep = "."), paste(ann$x, ann$y, sep = ".")); ann <- ann[ord, ]
pur <- sapply(sort(unique(base)), function(d) {
  a <- ann$annotation[base == d]; max(table(a)) / length(a) })

split_arb <- function(lab, d) { i <- which(lab == d)
  v <- if (diff(range(xy$x[i])) >= diff(range(xy$y[i]))) xy$x[i] else xy$y[i]
  lab[i[v > median(v)]] <- max(lab) + 1L; lab }
split_inf <- function(lab, d) { i <- which(lab == d); a <- ann$annotation[i]
  top <- names(sort(table(a), decreasing = TRUE))[1]
  lab[i[a != top]] <- max(lab) + 1L; lab }

d_arb <- order(-pur)[1]; d_inf <- order(pur)[1]
A <- split_arb(base, d_arb); I <- split_inf(base, d_inf)
da <- di <- numeric(0)
for (s in 1:8) {
  set.seed(2000 + s)
  fob <- setNames(sample(rep_len(1:NFOLD, length(ub))), ub); fold <<- fob[blk]
  b <- cv_deviance_explained(base)
  da <- c(da, cv_deviance_explained(A) - b)
  di <- c(di, cv_deviance_explained(I) - b)
}
cat(sprintf("\narbitrary split  (domain %d, purity %.2f): mean %+.5f  sd %.5f\n",
            d_arb, pur[d_arb], mean(da), sd(da)))
cat(sprintf("informative split(domain %d, purity %.2f): mean %+.5f  sd %.5f\n",
            d_inf, pur[d_inf], mean(di), sd(di)))
cat(sprintf("separation %.5f, or %.1f pooled sd\n", mean(di) - mean(da),
            (mean(di) - mean(da)) / sqrt((var(da) + var(di)) / 2)))
