## Diagnose stage-1 reproduction-check failures. For each gene: one local fit (published model),
## two independent draw sets (seeds A, B). Apply the same pass statistics to (A vs saved) and to
## (A vs B). If A-vs-B also fails, the failure is the check's assumption (Monte Carlo error
## independent across bins), not a refit discrepancy.
## Usage: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript diagnose-check.R Tlr1 Alox5 ...
suppressMessages({library(data.table); library(INLA); library(fmesher); library(splancs); library(inlabru); library(parallel)})
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
pdir <- file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects")
genes <- commandArgs(trailingOnly = TRUE); VB <- as.logical(Sys.getenv("VB", "TRUE")); inla.setOption(num.threads = "1:1")
sub <- readRDS(file.path(repo, "data-outputs/mouse-embryo-raw/controls/bin50-subregion-raw.rds"))[feat %in% genes]
sb <- as.points(matrix(c(15000, 25000, 25000, 15000, 15000, 11500, 11500, 18000, 18000, 11500), ncol = 2))
me <- diff(range(sb[, 1])) / 15; bo <- diff(range(sb[, 1])) / 3
mesh <- fm_mesh_2d(boundary = list(fm_segm(c(fm_segm(sb)))), max.edge = c(.15, 2) * me, offset = c(me, bo / 1.5), cutoff = me / 10, min.angle = 21)
matern <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, constr = TRUE, prior.range = c(500, .95), prior.sigma = c(.1, .05))
bru_options_set(control.compute = list(cpo = FALSE))
st <- function(la, lb, s, n = 1000) { d <- apply(la, 1, median) - apply(lb, 1, median); z <- d / (sqrt(2) * 1.2533 * s / sqrt(n))
  c(mean_z = mean(z), sd_z = sd(z), med_dlog = median(d), q99_abs = unname(quantile(abs(d), .99)),
    width_ratio = median((apply(la, 1, quantile, .975) - apply(la, 1, quantile, .025)) / (apply(lb, 1, quantile, .975) - apply(lb, 1, quantile, .025)))) }
res <- mclapply(genes, function(g) {
  rd <- sub[feat == g, .(x, y, feat.count, total.count)]; setkey(rd, x, y)
  fit <- bru(~ i(1) + f(cbind(x, y), model = matern),
             bru_obs(family = "poisson", data = rd, formula = feat.count ~ i + f, E = total.count), options = list(bru_verbose = 0, bru_max_iter = 1, control.inla = list(control.vb = list(enable = VB))))
  A <- log(generate(fit, rd, ~ exp(i + f), n.samples = 1000, seed = 11L)); B <- log(generate(fit, rd, ~ exp(i + f), n.samples = 1000, seed = 22L))
  sv <- readRDS(file.path(pdir, paste0(g, "-pred-poi.rds")))$lambda
  s <- (log(sv$q0.975) - log(sv$q0.025)) / (2 * qnorm(.975))
  S <- cbind(log(sv$q0.025), log(sv$median), log(sv$q0.975))
  ## saved summaries as a 3-column stand-in: only medians and widths are compared against saved
  dsv <- apply(A, 1, median) - log(sv$median)
  zsv <- dsv / (sqrt(2) * 1.2533 * s / sqrt(1000))
  wsv <- median((apply(A, 1, quantile, .975) - apply(A, 1, quantile, .025)) / (S[, 3] - S[, 1]))
  rbind(data.table(gene = g, comparison = "refit A vs saved", mean_z = mean(zsv), sd_z = sd(zsv), med_dlog = median(dsv),
                   q99_abs = quantile(abs(dsv), .99), width_ratio = wsv,
                   range = fit$summary.hyperpar[1, "mean"], sigma = fit$summary.hyperpar[2, "mean"]),
        data.table(gene = g, comparison = "refit A vs refit B", t(st(A, B, s)), range = NA, sigma = NA))
}, mc.cores = 3)
print(rbindlist(res), digits = 3)
