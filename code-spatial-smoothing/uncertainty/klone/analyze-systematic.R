## Phase 4 tabulation, run locally on the stage-2 outputs copied back from klone.
## Inputs: PAIRDIR (directory of <pair>.rds from uq-stage2.R), the published fields, the design.
## Outputs (in klone/results/): regional.csv (pair x 13 comparisons), crossfield.csv (field pairs
## sharing a gene x 10 regions), slices.csv (pair x 4 slices, summaries), summary.org.
## Run: PAIRDIR=<dir> OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript analyze-systematic.R
suppressMessages({library(data.table); library(excursions)})
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"; ud <- file.path(repo, "code-spatial-smoothing/uncertainty")
pdir <- Sys.getenv("PAIRDIR"); dir.create(file.path(ud, "klone/results"), showWarnings = FALSE)
TS <- c(1.25, 1.5, 2)
I <- readRDS(file.path(ud, "klone/stage/design-index.rds")); D <- readRDS(file.path(ud, "design/design-regions.rds"))
P <- fread(file.path(ud, "klone/stage/pairs.csv"))
mn <- readRDS(file.path(repo, "code-spatial-smoothing/controls/sweep/m-est-min-matrix.rds"))
keep <- !grepl("conflicted", mn$features); Mpub <- mn$counts[keep, ]; rownames(Mpub) <- mn$features[keep]
tc <- readRDS(file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
stopifnot(all(colnames(Mpub) == paste(tc$x, tc$y, sep = ".")))
Rpub <- sweep(Mpub, 2, tc$total.count, "/")
cmp <- rbind(data.table(r1 = c("box-heart-1", "box-cord-1", "box-heart-1"), r2 = c("box-heart-2", "box-cord-2", "box-liver")),
             as.data.table(t(combn(grep("^anat-", names(I$regions), value = TRUE), 2)))[, .(r1 = V1, r2 = V2)])
ci <- function(c) { q <- quantile(c, c(.025, .5, .975)); list(lo = exp(q[[1]]), med = exp(q[[2]]), hi = exp(q[[3]]), cred = q[[1]] > 0 || q[[3]] < 0) }
files <- file.path(pdir, paste0(P$pair, ".rds")); have <- file.exists(files)
cat("pairs with stage-2 output:", sum(have), "of", nrow(P), "\n")
X <- lapply(setNames(files[have], P$pair[have]), readRDS)
## soft flag (not an exclusion): pairs with a gene whose refit differs from the published fit
## beyond Monte Carlo noise (uq-check2.R); every summary is reported with and without them
ck <- fread(file.path(dirname(pdir), "check2.csv")); fg <- ck[flag == TRUE, gene]
P[, flagged := Lg %in% fg | Rg %in% fg]; fl <- P[flagged == TRUE, pair]
cat("flagged genes:", length(fg), " flagged pairs:", length(fl), "\n")

## ---- regional, same field
reg <- rbindlist(lapply(names(X), function(p) { x <- X[[p]]; pl <- Rpub[p, ]
  rbindlist(lapply(seq_len(nrow(cmp)), function(i) { r1 <- cmp$r1[i]; r2 <- cmp$r2[i]; k1 <- I$regions[[r1]]; k2 <- I$regions[[r2]]
    a <- mean(pl[k1]); b <- mean(pl[k2])
    j <- ci(log(x$region_mean[r1, ]) - log(x$region_mean[r2, ])); d <- ci(log(x$region_mean_decor[r1, ]) - log(x$region_mean_decor[r2, ]))
    data.table(pair = p, r1, r2, plug_ratio = a / b, wilcox_p = suppressWarnings(wilcox.test(pl[k1], pl[k2])$p.value),
               post_med = j$med, lo = j$lo, hi = j$hi, credible = j$cred, dec_lo = d$lo, dec_hi = d$hi, dec_credible = d$cred) })) }))
reg[, flagged := pair %in% fl]
fwrite(reg, file.path(ud, "klone/results/regional.csv"))

## ---- regional, cross field (field pairs sharing a ligand or a receptor), shape-normalized
lig <- P[have, .(pair, g = Lg)]; rec <- P[have, .(pair, g = Rg)]
fp <- unique(rbind(lig[lig, on = "g", allow.cartesian = TRUE, .(a = pair, b = i.pair, shared = "ligand")],
                   rec[rec, on = "g", allow.cartesian = TRUE, .(a = pair, b = i.pair, shared = "receptor")])[a < b], by = c("a", "b"))
cat("cross-field pairs:", nrow(fp), "\n")
regs <- names(I$regions)
shp <- lapply(X, function(x) log(sweep(x$region_mean, 2, x$tissue_mean, "/")))             # 10 x 1000
pls <- lapply(setNames(names(X), names(X)), function(p) Rpub[p, ] / mean(Rpub[p, ]))
plr <- sapply(names(X), function(p) sapply(regs, function(r) mean(pls[[p]][I$regions[[r]]])))   # regions x pairs
crx <- rbindlist(lapply(seq_len(nrow(fp)), function(i) { a <- fp$a[i]; b <- fp$b[i]; da <- shp[[a]] - shp[[b]]
  rbindlist(lapply(regs, function(r) { k <- I$regions[[r]]; j <- ci(da[r, ])
    data.table(a, b, shared = fp$shared[i], region = r, plug_ratio = plr[r, a] / plr[r, b],
               wilcox_signed_p = suppressWarnings(wilcox.test(pls[[a]][k], pls[[b]][k], paired = TRUE)$p.value),
               post_med = j$med, lo = j$lo, hi = j$hi, credible = j$cred) })) }))
crx[, flagged := a %in% fl | b %in% fl]
fwrite(crx, file.path(ud, "klone/results/crossfield.csv"))

## ---- slices, same field: deviation from the slice average; joint, simultaneous, pointwise; labels
off <- cumsum(c(0, lengths(I$slices)))
sls <- rbindlist(lapply(names(X), function(p) { x <- X[[p]]; L <- x$slice_logM_i16 / 1000 - 20
  rbindlist(lapply(seq_along(I$slices), function(s) { rr <- (off[s] + 1):off[s + 1]; Lm <- exp(L[rr, , drop = FALSE])
    dv <- log(Lm) - matrix(log(colMeans(Lm)), length(rr), ncol(Lm), byrow = TRUE)
    pe <- excursions.mc(dv, alpha = .05, u = 0, type = ">"); ne <- excursions.mc(dv, alpha = .05, u = 0, type = "<")
    invisible(capture.output(sc <- simconf.mc(dv, alpha = .05)))
    pl <- Rpub[p, I$slices[[s]]]; pdv <- log(pl) - log(mean(pl)); inj <- as.logical(pe$E) | as.logical(ne$E)
    pw <- rowMeans(dv > 0) >= .95 | rowMeans(dv < 0) >= .95
    data.table(pair = p, slice = names(I$slices)[s], n = length(rr), n_pointwise = sum(pw), n_simband = sum(sc$a > 0 | sc$b < 0),
               n_joint = sum(inj), ld_125 = sum(abs(pdv) >= log(1.25)), ld_125_joint = sum(abs(pdv) >= log(1.25) & inj),
               ld_15 = sum(abs(pdv) >= log(1.5)), ld_15_joint = sum(abs(pdv) >= log(1.5) & inj),
               ld_2 = sum(abs(pdv) >= log(2)), ld_2_joint = sum(abs(pdv) >= log(2) & inj),
               alike_joint_125 = sum(abs(pdv) < log(1.25) & inj)) })) }))
sls[, flagged := pair %in% fl]
fwrite(sls, file.path(ud, "klone/results/slices.csv"))

## ---- summary
tabf <- function(d, T0) d[, .(looks_diff = sum(pmax(plug_ratio, 1 / plug_ratio) >= T0),
  looks_diff_not_credible = sum(pmax(plug_ratio, 1 / plug_ratio) >= T0 & !credible),
  looks_alike_credible = sum(pmax(plug_ratio, 1 / plug_ratio) < T0 & credible), credible = sum(credible), n = .N)]
sink(file.path(ud, "klone/results/summary.org"))
cat("#+TITLE: Systematic run summary (generated by analyze-systematic.R)\n\n")
for (T0 in TS) { cat(sprintf("* T = %.2f\n\nRegional, same field (all pairs; then excluding flagged pairs):\n", T0)); print(tabf(reg, T0)); print(tabf(reg[flagged == FALSE], T0))
  cat("\nWilcoxon p < 0.05:", reg[, sum(wilcox_p < .05)], "  decorrelated credible:", reg[, sum(dec_credible)], "\n\nCross field (all; then excluding flagged):\n"); print(tabf(crx, T0)); print(tabf(crx[flagged == FALSE], T0)); cat("\n") }
cat("* Slices (summed over pairs)\n\n"); print(sls[, lapply(.SD, sum), by = slice, .SDcols = n:alike_joint_125])
sink()
cat(readLines(file.path(ud, "klone/results/summary.org")), sep = "\n")
