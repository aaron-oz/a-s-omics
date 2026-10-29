## Prior-sensitivity refits (REVIEW-2026-10-29.org item P2; logged in DESIGN.org 2026-10-29).
## Identical to refit-draws.R except the penalized-complexity prior, chosen by env PRIOR:
##   a: P(range < 250) = 0.95 (range prior halved), sigma prior unchanged
##   b: P(range < 1000) = 0.95 (range prior doubled), sigma prior unchanged
##   c: P(sigma > 1) = 0.05 (sigma prior relaxed), range prior unchanged
## Output: /var/home/aoz/data/a-s-omics-uncertainty/draws-prior-<PRIOR>/<Gene>-lambda-draws.rds
## Usage: PRIOR=a NCORE=6 OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript refit-draws-prior.R Calr Itga2b ...
suppressMessages({library(data.table); library(INLA); library(fmesher); library(splancs)
                  library(inlabru); library(parallel)})
repo  <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
PRIOR <- Sys.getenv("PRIOR"); stopifnot(PRIOR %in% c("a", "b", "c"))
odir  <- paste0("/var/home/aoz/data/a-s-omics-uncertainty/draws-prior-", PRIOR); dir.create(odir, showWarnings = FALSE)
genes <- commandArgs(trailingOnly = TRUE)
ncore <- as.integer(Sys.getenv("NCORE", "1")); NSAMP <- as.integer(Sys.getenv("NSAMP", "1000"))
inla.setOption(num.threads = "1:1")
sub <- readRDS(file.path(repo, "data-outputs/mouse-embryo-raw/controls/bin50-subregion-raw.rds"))
sub <- sub[feat %in% genes]; stopifnot(all(genes %in% sub$feat))
sb <- as.points(matrix(c(15000, 25000, 25000, 15000, 15000, 11500, 11500, 18000, 18000, 11500), ncol = 2))
max.edge <- diff(range(sb[, 1])) / (3 * 5); bound.outer <- diff(range(sb[, 1])) / 3
mesh <- fm_mesh_2d(boundary = list(fm_segm(c(fm_segm(sb)))), max.edge = c(.15, 2) * max.edge,
                   offset = c(max.edge, bound.outer / 1.5), cutoff = max.edge / 10, min.angle = 21)
pr <- switch(PRIOR, a = list(range = c(250, .95), sigma = c(.1, .05)), b = list(range = c(1000, .95), sigma = c(.1, .05)),
             c = list(range = c(500, .95), sigma = c(1, .05)))
matern <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, constr = TRUE, prior.range = pr$range, prior.sigma = pr$sigma)
bru_options_set(control.compute = list(cpo = FALSE))
fit_one <- function(g) {
  t0 <- Sys.time()
  rd <- sub[feat == g, .(x, y, feat.count, total.count)]; setkey(rd, x, y)
  comps <- ~ feat.count.int(1) + feat.count.field(cbind(x, y), model = matern)
  lik <- bru_obs(family = "poisson", data = rd, formula = feat.count ~ feat.count.int + feat.count.field, E = total.count)
  fit <- bru(comps, lik, options = list(bru_verbose = 0, bru_max_iter = 1))
  seed <- 20261029L + match(g, genes)
  dr <- generate(fit, rd, ~ exp(feat.count.int + feat.count.field), n.samples = NSAMP, seed = seed)
  out <- list(gene = g, prior = PRIOR, pc = pr, xy = as.matrix(rd[, .(x, y)]), draws = dr, seed = seed,
              hyper = fit$summary.hyperpar, seconds = as.numeric(Sys.time() - t0, units = "secs"))
  saveRDS(out, file.path(odir, paste0(g, "-lambda-draws.rds")))
  sprintf("%s [%s]: range %.0f sigma %.3f, %.0f s", g, PRIOR, fit$summary.hyperpar[1, "mean"], fit$summary.hyperpar[2, "mean"], out$seconds)
}
res <- mclapply(genes, function(g) tryCatch(fit_one(g), error = function(e) paste(g, "FAILED:", conditionMessage(e))), mc.cores = ncore)
cat(unlist(res), sep = "\n")
