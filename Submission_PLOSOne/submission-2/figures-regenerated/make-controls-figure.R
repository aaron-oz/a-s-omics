## Build the controls figure (S10 Fig) for reviewer comment R2.8.
##
## Four panels:
##   A  spatial permutation test: null distribution of ELSA(d50) over 49 coordinate
##      permutations, with the observed value marked
##   B  the same for the number of domains recovered
##   C  emergent signalling domains beside the published anatomical annotation, same section
##   D  concordance against a 499-permutation null: adjusted Rand index, observed marked
##
## Run: distrobox enter emacs-r -- Rscript make-controls-figure.R
suppressMessages({library(data.table); library(fields)})
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
cd   <- file.path(repo, "data-outputs/mouse-embryo-raw/controls")
odir <- file.path(repo, "Submission_PLOSOne/submission-2/figures-regenerated")

B  <- rbindlist(lapply(Sys.glob(file.path(cd, "control-B-klone", "control-B-perm*.csv")), fread))
obs <- B[perm_id == 0]; nul <- B[perm_id != 0]
C  <- readRDS(file.path(cd, "control-C-result.rds"))
lab <- readRDS(file.path(cd, "observed-domain-labels.rds"))
ann <- fread(file.path(cd, "e16-fov-annotation.csv"))

png(file.path(odir, "S10_Fig_controls.png"), width = 13, height = 9.5, units = "in", res = 300)
layout(matrix(c(1,2, 3,4, 5,5), nrow = 3, byrow = TRUE), heights = c(1, 1.3, 1))
par(mar = c(4.2, 4.4, 3.2, 1.4), cex.lab = 1.05)

## A ---- ELSA null
h <- hist(nul$elsa_d50, breaks = 12, plot = FALSE)
plot(h, col = "grey80", border = "white", xlim = range(c(obs$elsa_d50, nul$elsa_d50)) * c(0.85, 1.05),
     main = "A  Spatial permutation test: domain coherence", xlab = "mean ELSA at radius 50", ylab = "permutations")
abline(v = obs$elsa_d50, col = "firebrick", lwd = 3)
text(obs$elsa_d50, max(h$counts) * 0.92, sprintf("  observed %.3f", obs$elsa_d50),
     col = "firebrick", adj = 0, font = 2)
text(mean(nul$elsa_d50), max(h$counts) * 0.55,
     sprintf("49 permutations\nmean %.3f", mean(nul$elsa_d50)), adj = 0.5, col = "grey25")
mtext("lower = more spatially coherent", side = 3, line = 0.1, cex = 0.72, col = "grey35")

## B ---- domain count null
h2 <- hist(nul$n_clusters, breaks = seq(min(nul$n_clusters) - 0.5, max(nul$n_clusters) + 0.5, 1), plot = FALSE)
plot(h2, col = "grey80", border = "white", xlim = c(min(nul$n_clusters) - 1, obs$n_clusters + 2),
     main = "B  Spatial permutation test: number of domains", xlab = "domains recovered", ylab = "permutations")
abline(v = obs$n_clusters, col = "firebrick", lwd = 3)
text(obs$n_clusters, max(h2$counts) * 0.92, sprintf("  observed %d", obs$n_clusters),
     col = "firebrick", adj = 0, font = 2)

## C ---- domains vs anatomy, side by side
setnames(lab, "domain", "dom")
D <- merge(ann, lab, by = c("x","y"))
nx <- length(unique(D$x)); ny <- length(unique(D$y))
par(mar = c(3.4, 3.4, 3.2, 1.0))
quilt.plot(D$x, D$y, D$dom, nx = nx, ny = ny, asp = 1, add.legend = FALSE,
           col = hcl.colors(length(unique(D$dom)), "Dark 3"),
           main = "C  Emergent signalling domains", xlab = "", ylab = "")
al <- as.integer(factor(D$annotation))
quilt.plot(D$x, D$y, al, nx = nx, ny = ny, asp = 1, add.legend = FALSE,
           col = hcl.colors(length(unique(al)), "Set 3"),
           main = "      Published anatomical annotation, same section", xlab = "", ylab = "")

## D ---- ARI null
par(mar = c(4.2, 4.4, 3.2, 1.4))
## bins span the whole plotted range, so the null is a visible bar rather than a sub-pixel spike
brk <- seq(min(C$null) - 0.005, C$ari * 1.06, length.out = 45)
h3 <- hist(C$null, breaks = brk, plot = FALSE)
plot(h3, col = "grey80", border = "grey60",
     main = "D  Concordance with the anatomical annotation",
     xlab = "adjusted Rand index", ylab = "permutations")
abline(v = C$ari, col = "firebrick", lwd = 3)
text(C$ari, max(h3$counts) * 0.88, sprintf("observed %.3f  ", C$ari), col = "firebrick", adj = 1, font = 2)
text(max(C$null) + 0.02, max(h3$counts) * 0.55,
     sprintf("499 permutations,\nall below %.4f", max(C$null)), adj = 0, col = "grey25")
mtext(sprintf("p = %.3f", C$p), side = 3, line = 0.1, cex = 0.72, col = "grey35")
dev.off()
cat("wrote S10_Fig_controls.png\n")
cat(sprintf("  A/B: observed ELSA %.3f, %d domains; null mean %.3f, %d-%d domains\n",
            obs$elsa_d50, obs$n_clusters, mean(nul$elsa_d50), min(nul$n_clusters), max(nul$n_clusters)))
cat(sprintf("  D:   observed ARI %.3f vs null max %.4f\n", C$ari, max(C$null)))
