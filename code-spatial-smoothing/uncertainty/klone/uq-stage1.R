## Systematic run, stage 1 (one array task = one batch of genes). Same model as pilot/refit-draws.R
## (see pilot/FIDELITY.org). Per gene: fit, 1,000 joint draws of lambda at all 26,000 bins, an
## on-cluster reproduction check against the saved medians (2,000 fixed bins), and the draws
## written as int16: q = round((log lambda + 20) * 1000), i.e. log-scale resolution 0.001.
## Idempotent: a gene whose .i16 and check row already exist is skipped; files are written to a
## temporary name and renamed, so a preempted task leaves no partial file.
## env: ROOT (input dir), DRAWS (output dir for draws), OUT (dir for small csv), BATCH (genes per task)
suppressMessages({library(data.table); library(INLA); library(fmesher); library(splancs); library(inlabru)})
root <- Sys.getenv("ROOT"); ddir <- Sys.getenv("DRAWS"); odir <- Sys.getenv("OUT")
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0")); B <- as.integer(Sys.getenv("BATCH", "10"))
NSAMP <- 1000L
inla.setOption(num.threads = "1:1")
genes <- readLines(file.path(root, "genes.txt"))
mine <- genes[(task * B + 1):min(length(genes), (task + 1) * B)]
mine <- mine[!is.na(mine)]
sub <- readRDS(file.path(root, "bin50-subregion-raw.rds")); sub <- sub[feat %in% mine]
chk <- readRDS(file.path(root, "check-bins.rds"))
sb <- as.points(matrix(c(15000, 25000, 25000, 15000, 15000, 11500, 11500, 18000, 18000, 11500), ncol = 2))
max.edge <- diff(range(sb[, 1])) / (3 * 5); bound.outer <- diff(range(sb[, 1])) / 3
mesh <- fm_mesh_2d(boundary = list(fm_segm(c(fm_segm(sb)))), max.edge = c(.15, 2) * max.edge,
                   offset = c(max.edge, bound.outer / 1.5), cutoff = max.edge / 10, min.angle = 21)
matern <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, constr = TRUE, prior.range = c(500, .95), prior.sigma = c(.1, .05))
bru_options_set(control.compute = list(cpo = FALSE))
for (g in mine) {
  f16 <- file.path(ddir, paste0(g, ".i16")); fck <- file.path(odir, "check", paste0(g, ".csv"))
  if (file.exists(f16) && file.exists(fck)) next
  t0 <- Sys.time()
  rd <- sub[feat == g, .(x, y, feat.count, total.count)]; setkey(rd, x, y); stopifnot(nrow(rd) == 26000)
  comps <- ~ feat.count.int(1) + feat.count.field(cbind(x, y), model = matern)
  lik <- bru_obs(family = "poisson", data = rd, formula = feat.count ~ feat.count.int + feat.count.field, E = total.count)
  fit <- bru(comps, lik, options = list(bru_verbose = 0, bru_max_iter = 1))
  seed <- 20261021L + match(g, genes)
  lg <- log(generate(fit, rd, ~ exp(feat.count.int + feat.count.field), n.samples = NSAMP, seed = seed))
  q <- round((lg + 20) * 1000); stopifnot(all(q >= -32768 & q <= 32767))
  ## reproduction check on the fixed bins, from the quantized draws (so it also checks quantization)
  lq <- q[chk$bins, ] / 1000 - 20; sv <- chk$saved[[g]]
  s <- (log(sv[, 3]) - log(sv[, 2])) / (2 * qnorm(.975))
  dlog <- apply(lq, 1, median) - log(sv[, 1]); z <- dlog / (sqrt(2) * 1.2533 * s / sqrt(NSAMP))
  row <- data.table(gene = g, seed = seed, mean_z = mean(z), sd_z = sd(z), med_dlog = median(dlog),
                    pass = abs(mean(z)) < 0.2 && sd(z) < 1.3 && abs(median(dlog)) < 0.02,
                    quant_err_max = max(abs(lg[chk$bins, ] - lq)), range = fit$summary.hyperpar[1, "mean"],
                    sigma = fit$summary.hyperpar[2, "mean"], secs = round(as.numeric(Sys.time() - t0, units = "secs")),
                    inla = as.character(packageVersion("INLA")), inlabru = as.character(packageVersion("inlabru")))
  tmp <- paste0(f16, ".tmp"); con <- file(tmp, "wb"); writeBin(as.integer(q), con, size = 2); close(con)
  file.rename(tmp, f16)
  dir.create(dirname(fck), showWarnings = FALSE, recursive = TRUE)
  fwrite(row, paste0(fck, ".tmp")); file.rename(paste0(fck, ".tmp"), fck)
  message(sprintf("%s: %s s, sd_z %.3f, pass %s", g, row$secs, row$sd_z, row$pass))
}
