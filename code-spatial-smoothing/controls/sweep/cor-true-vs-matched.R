## Are true ligand-receptor partners more spatially correlated than matched random ones?
##
## AOZ's hypothesis for the matched-dominance scramble result: true signalling partners
## should have correlated spatial fields, and that alone may explain both "wrong
## direction" findings. For jointly normal (L,R) with correlation rho, writing
## U=(L+R)/2 and V=(L-R)/2, min = U - |V|; U and V are independent and |V| is
## uncorrelated with every linear function of (L,R), so the best linear predictor of the
## minimum is exactly U and
##      R2(rho) = (1+rho) / [ (1+rho) + (1-rho)(1-2/pi) ]
## which is monotone INCREASING in rho. Higher partner correlation => less nonlinear
## content, and (because the structured parts cancel in the difference) noisier switching.
## Both observed effects follow if rho_true > rho_null.
##
## Uses the exact matched partners chosen by matched-dominance-scramble.R, so the
## comparison is paired on the same mechanisms.

gf <- readRDS("gene-field-matrix.rds"); G <- gf$G
ms <- read.csv("matched-dominance-scramble.csv", stringsAsFactors = FALSE)
li <- match(ms$ligand, colnames(G)); ri <- match(ms$receptor, colnames(G))
mi <- match(ms$matched_partner, colnames(G))
stopifnot(!any(is.na(li)), !any(is.na(ri)), !any(is.na(mi)))
n <- nrow(ms); cat("paired mechanisms: ", n, "\n\n", sep = "")

rho_t <- rho_n <- sp_t <- sp_n <- numeric(n)
for (i in seq_len(n)) {
  L <- G[, li[i]]
  rho_t[i] <- cor(L, G[, ri[i]]);            rho_n[i] <- cor(L, G[, mi[i]])
  sp_t[i]  <- cor(L, G[, ri[i]], method="spearman")
  sp_n[i]  <- cor(L, G[, mi[i]], method="spearman")
}

f <- function(x) sprintf("%.4f", quantile(x, c(.1,.25,.5,.75,.9), na.rm=TRUE))
line <- function(lab, x) cat(sprintf("  %-20s mean %+.4f  median %+.4f   deciles %s\n",
  lab, mean(x), median(x), paste(f(x), collapse=" ")))

cat("PEARSON cor(L, partner)\n"); line("true partner", rho_t); line("matched random", rho_n)
d <- rho_t - rho_n
cat(sprintf("  paired difference: mean %+.4f  median %+.4f  sd %.4f\n", mean(d), median(d), sd(d)))
cat(sprintf("  true higher in %.1f%% of mechanisms; paired t p = %.3g; Wilcoxon p = %.3g\n\n",
  100*mean(d > 0), t.test(rho_t, rho_n, paired=TRUE)$p.value,
  wilcox.test(rho_t, rho_n, paired=TRUE)$p.value))

cat("SPEARMAN cor(L, partner)\n"); line("true partner", sp_t); line("matched random", sp_n)
d2 <- sp_t - sp_n
cat(sprintf("  paired difference: mean %+.4f  median %+.4f\n", mean(d2), median(d2)))
cat(sprintf("  true higher in %.1f%%; Wilcoxon p = %.3g\n\n", 100*mean(d2 > 0),
  wilcox.test(sp_t, sp_n, paired=TRUE)$p.value))

## does the correlation gap ACCOUNT for the linear-R2 gap? (mediation check)
R2t <- ms$R2_true; R2n <- ms$R2_null
cat("MEDIATION: does rho explain the linear-R^2 difference?\n")
cat(sprintf("  raw R^2 gap (true - null): mean %+.5f\n", mean(R2t - R2n, na.rm=TRUE)))
long <- data.frame(R2 = c(R2t, R2n), rho = c(rho_t, rho_n),
                   dom = c(ms$dom_true, ms$dom_null),
                   arm = rep(c(1, 0), each = n))
m0 <- lm(R2 ~ arm, long); m1 <- lm(R2 ~ arm + rho + dom, long)
cat(sprintf("  arm coefficient, unadjusted:            %+.5f (p = %.3g)\n",
  coef(m0)["arm"], summary(m0)$coefficients["arm", 4]))
cat(sprintf("  arm coefficient, adjusted for rho+dom:  %+.5f (p = %.3g)\n",
  coef(m1)["arm"], summary(m1)$coefficients["arm", 4]))
cat(sprintf("  attenuation: %.1f%% of the arm effect is accounted for by rho and dominance\n",
  100 * (1 - coef(m1)["arm"] / coef(m0)["arm"])))

## sanity: is the empirical R2-vs-rho relationship increasing, as the theory says?
cat(sprintf("\n  cor(rho, linear R^2) across all %d observations: %+.3f (theory: positive)\n",
  nrow(long), cor(long$rho, long$R2, use="complete.obs")))
pred <- function(r) (1+r) / ((1+r) + (1-r)*(1-2/pi))
cat(sprintf("  bivariate-normal predicted R^2 at mean rho_true %.4f: %.4f\n", mean(rho_t), pred(mean(rho_t))))
cat(sprintf("  bivariate-normal predicted R^2 at mean rho_null %.4f: %.4f\n", mean(rho_n), pred(mean(rho_n))))
cat(sprintf("  predicted gap %+.5f  vs observed gap %+.5f\n",
  pred(mean(rho_t)) - pred(mean(rho_n)), mean(R2t - R2n, na.rm=TRUE)))

write.csv(data.frame(mechanism=ms$mechanism, rho_true=rho_t, rho_null=rho_n,
  spearman_true=sp_t, spearman_null=sp_n), "cor-true-vs-matched.csv", row.names=FALSE)
cat("\nwrote cor-true-vs-matched.csv\n")
