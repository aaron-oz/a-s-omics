## Phase 3 pilot: refit genes with the published Poisson model and save JOINT posterior draws
## of the rate field lambda(s) at all 26,000 bins. Model specification copied from
## raw-counts/02-array-run-model.R (Poisson branch) and raw-counts/01-array-run-job.R (mesh);
## see FIDELITY.org for the row-by-row comparison.
## Draws: inlabru generate(), D = 1000, fixed seed per gene, INLA single-threaded.
## Output (outside Dropbox): /var/home/aoz/data/a-s-omics-uncertainty/draws/<Gene>-lambda-draws.rds
##   list(gene, xy (26000 x 2, sorted by x then y), draws (26000 x 1000 matrix of lambda), seed,
##        hyper (posterior summary of hyperparameters), versions, seconds)
## Usage: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript refit-draws.R Cd44 [Wnt4 ...]
## Env: NCORE (parallel genes, default 1), NSAMP (default 1000)
suppressMessages({library(data.table); library(INLA); library(fmesher); library(splancs)
                  library(inlabru); library(parallel)})
repo  <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
odir  <- "/var/home/aoz/data/a-s-omics-uncertainty/draws"
genes <- commandArgs(trailingOnly = TRUE)
ncore <- as.integer(Sys.getenv("NCORE", "1")); NSAMP <- as.integer(Sys.getenv("NSAMP", "1000"))
inla.setOption(num.threads = "1:1")

sub <- readRDS(file.path(repo, "data-outputs/mouse-embryo-raw/controls/bin50-subregion-raw.rds"))
sub <- sub[feat %in% genes]
stopifnot(all(genes %in% sub$feat))

## mesh: identical construction to 01-array-run-job.R with agg.size = 50 (inner multiplier 0.15)
sb <- as.points(matrix(c(15000, 25000, 25000, 15000, 15000, 11500, 11500, 18000, 18000, 11500), ncol = 2))
max.edge <- diff(range(sb[, 1])) / (3 * 5); bound.outer <- diff(range(sb[, 1])) / 3
mesh <- fm_mesh_2d(boundary = list(fm_segm(c(fm_segm(sb)))), max.edge = c(.15, 2) * max.edge,
                   offset = c(max.edge, bound.outer / 1.5), cutoff = max.edge / 10, min.angle = 21)
## prior: the Poisson branch of 02-array-run-model.R (lines 63-67)
matern <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, constr = TRUE,
                              prior.range = c(500, .95), prior.sigma = c(.1, .05))
bru_options_set(control.compute = list(cpo = FALSE))

fit_one <- function(g) {
  t0 <- Sys.time()
  rd <- sub[feat == g, .(x, y, feat.count, total.count)]; setkey(rd, x, y)
  comps <- ~ feat.count.int(1) + feat.count.field(cbind(x, y), model = matern)
  lik <- bru_obs(family = "poisson", data = rd,
                 formula = feat.count ~ feat.count.int + feat.count.field, E = total.count)
  fit <- bru(comps, lik, options = list(bru_verbose = 0, bru_max_iter = 1))
  seed <- 20261016L + match(g, genes)
  dr <- generate(fit, rd, ~ exp(feat.count.int + feat.count.field), n.samples = NSAMP, seed = seed)
  out <- list(gene = g, xy = as.matrix(rd[, .(x, y)]), draws = dr, seed = seed,
              hyper = fit$summary.hyperpar,
              versions = c(INLA = as.character(packageVersion("INLA")),
                           inlabru = as.character(packageVersion("inlabru")),
                           fmesher = as.character(packageVersion("fmesher"))),
              seconds = as.numeric(Sys.time() - t0, units = "secs"))
  saveRDS(out, file.path(odir, paste0(g, "-lambda-draws.rds")))
  sprintf("%s: %d bins x %d draws, %.0f s", g, nrow(dr), ncol(dr), out$seconds)
}
res <- mclapply(genes, function(g) tryCatch(fit_one(g), error = function(e) paste(g, "FAILED:", conditionMessage(e))),
                mc.cores = ncore)
cat(unlist(res), sep = "\n")
