## Noise floor: how much does the score move from the fold assignment alone?
## Source the core ONCE (it sets its own seed), then redraw folds in place.
S <- "/tmp/claude-1000/-var-home-aoz-Dropbox-genetics-a-s-omics/49644765-1ca5-4f64-867d-f5bd061fc895/scratchpad/cvcore.R"
source(S)
ROOT <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
e <- new.env()
load(file.path(ROOT, "sam_sandbox/clustering vs features/m.est.min",
               "clustering vs features experiment",
               "meta.data.feature.depth.experiment.2025-07-12.Robj"), envir = e)
md <- get(ls(e)[1], envir = e)
cc <- grep("featureClusters$", colnames(md), value = TRUE)
base <- as.integer(as.factor(md[[cc[which(as.integer(sub("featureClusters$","",cc)) == 1481)]]]))

sc <- numeric(0)
for (s in 1:8) {
  set.seed(1000 + s)
  fob <- setNames(sample(rep_len(1:NFOLD, length(ub))), ub)
  fold <<- fob[blk]
  sc <- c(sc, cv_deviance_explained(base))
  message(sprintf("  fold seed %d: %.5f", s, tail(sc, 1)))
}
cat(sprintf("\nbaseline over 8 fold assignments: mean %.5f  sd %.5f  range %.5f\n",
            mean(sc), sd(sc), diff(range(sc))))
cat("calibration deltas for comparison: arbitrary splits -0.00023 to -0.00039;",
    "informative splits -0.00006 to +0.00194\n")
