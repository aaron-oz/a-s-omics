## Fig 2, data parameters, redrawn in code (2026-10-23) from the original hand-made image
## (plos-port/img/Data Parameters.png): same panels, white background, no text below 8 pt.
## Ligand and receptor sets: FANTOM5 ligands and receptors among the 912 modeled genes (a gene
## that is both counts in both). Whole-section totals: bin-50 sums of the raw GEM file
## (data-inputs/mouse-embryo-raw/full-raw-t16-data.tsv), made by the awk one-liner recorded below.
##   awk -F'\t' 'NR>7 { k = int($2/50) "\t" int($3/50); s[k] += $4 } END { for (k in s) print k "\t" s[k] }' \
##     full-raw-t16-data.tsv > whole-section-bin50-totals.tsv
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript fig2-data-parameters.R
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(ragg)})
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"; FONT <- "Liberation Sans"; PT <- 8
sub <- readRDS(file.path(repo, "data-outputs/mouse-embryo-raw/controls/bin50-subregion-raw.rds"))
e <- new.env(); load(file.path(repo, "code-spatial-smoothing/controls/sweep/fantom-hierarchy.Robj"), envir = e)
fan <- as.data.table(get("output", e))
G <- unique(sub$feat); up <- toupper(G)
lig <- G[up %in% fan$LIGAND]; rec <- G[up %in% fan$RECEPTOR]
cat("modeled genes:", length(G), " ligands:", length(lig), " receptors:", length(rec), "\n")
tot <- unique(sub[, .(x, y, total.count)])
agg <- function(gs) sub[feat %in% gs, .(cnt = sum(feat.count), nfeat = sum(feat.count > 0)), by = .(x, y)]
L <- agg(lig); R <- agg(rec)
d <- merge(merge(tot, L, by = c("x", "y")), R, by = c("x", "y"), suffixes = c(".L", ".R"))
d[, `:=`(i = (x - 15025) / 50 + 1, j = (y - 11525) / 50 + 1, fL = cnt.L / total.count, fR = cnt.R / total.count)]
W <- fread("/var/home/aoz/data/a-s-omics-uncertainty/fig2/whole-section-bin50-totals.tsv", col.names = c("bx", "by", "n"))
## display mask only: bins under 1,000 counts are bare chip (log10 histogram is bimodal, modes near 300 and 5,000)
W <- W[n >= 1000]
th <- theme_void(base_size = PT, base_family = FONT) +
  theme(plot.title = element_text(size = PT + 1, hjust = 0), legend.title = element_blank(), legend.text = element_text(size = PT),
        legend.position = "bottom", legend.key.width = unit(0.45, "in"), legend.key.height = unit(0.07, "in"),
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(3, 4, 3, 4))
mp <- function(df, v, t, pal, lim = NULL) ggplot(df, aes(i, j, fill = .data[[v]])) + geom_raster() + coord_fixed(expand = FALSE) +
  scale_fill_viridis_c(option = pal, limits = lim, oob = scales::squish, n.breaks = 4) + th + labs(title = t)
q <- function(v) quantile(v, c(0.005, 0.995))
pW <- ggplot(W, aes(bx, by, fill = pmin(n, quantile(n, 0.995)))) + geom_raster() + coord_fixed(expand = FALSE) +
  scale_fill_viridis_c(option = "cividis") +
  annotate("rect", xmin = 300, xmax = 500, ymin = 230, ymax = 360, fill = NA, colour = "white", linetype = "22", linewidth = 0.5) +
  th + labs(title = "A  Total counts, whole section")
pF <- mp(d, "total.count", "B  Total counts, field of view", "cividis", q(d$total.count))
lc <- range(q(c(d$cnt.L, d$cnt.R))); lf <- range(q(c(d$nfeat.L, d$nfeat.R))); lr <- range(q(c(d$fL, d$fR)))
right <- (mp(d, "cnt.L", "C  Ligand counts", "inferno", lc) | mp(d, "cnt.R", "D  Receptor counts", "inferno", lc)) /
  (mp(d, "nfeat.L", "E  Ligand genes detected", "viridis", lf) | mp(d, "nfeat.R", "F  Receptor genes detected", "viridis", lf)) /
  (mp(d, "fL", "G  Ligand share of counts", "mako", lr) | mp(d, "fR", "H  Receptor share of counts", "mako", lr))
g <- (pW / pF) | right
agg_png("fig2-data-parameters.png", width = 7.5, height = 5.6, units = "in", res = 600); print(g + plot_layout(widths = c(1, 2))); invisible(dev.off())
cat("ok\n")
