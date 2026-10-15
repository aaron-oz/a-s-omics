## Rebuild Figure 6 (feature-depth experiment) at the PLOS page size, correct aspect, and
## legible type.
##
## Why: the published Fig 6 has the same defect R2.6 flagged in Figs 3 and 4 (its spatial
## panels come from ggplot in Sup7 without coord_fixed(), so they are stretched to the device),
## its baked-in labels counted 1481 columns, and at the PLOS maximum height its text would
## print at roughly 3 to 4 pt against the 8 pt minimum. It also had no layout script.
##
## Nothing is re-clustered. Every panel is drawn from the stored outputs of the 2025-07-12
## m.est.min run (cluster labels, UMAP embeddings, per-bin ELSA, summary trends), and the
## random-subset curve from the stored Seurat replication results.
##
## Layout, transposed like Fig 5: depths as columns (10, 150, 650, 1300, all 1,477 pairs);
## rows are UMAP, domains in tissue space, entrogram; then one row of summary curves with the
## five depths marked; then the eight 3D surfaces, reused as an image from the published
## figure (they are perspective renders, unaffected by the aspect defect).
##
## Run: distrobox enter emacs-r -- Rscript regenerate-fig6.R
## Output: Fig6-composite.png, 7.5 in wide, 600 dpi.

suppressMessages({
  library(Seurat); library(ggplot2); library(data.table); library(patchwork); library(png)
})

repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
src  <- file.path(repo, "sam_sandbox/clustering vs features/m.est.min")
sw   <- file.path(repo, "code-spatial-smoothing/controls/sweep")
odir <- file.path(repo, "Submission_PLOSOne/submission-2/figures-regenerated")
pub_img <- file.path(odir, "Fig6-relabeled.png")   # the published figure, labels corrected

DEPTHS <- c(10, 150, 650, 1300, 1481)
lab_depth <- function(n) if (n == 1481) "all 1,477 pairs" else sprintf("%s pairs", format(n, big.mark = ","))
FEAT_NUM <- c(10, 50, seq(100, 1450, 50), 1481)    # the published sweep, in stored order

## published palette (Sup6 / Sup7); cluster k is drawn in colour k + 1, as in Sup7
cols.global <- unique(c(
  'firebrick','steelblue','springgreen','purple','salmon','skyblue','navyblue','violetred',
  '#40376E','violetred','grey20','tomato','sandybrown','saddlebrown',
  'royalblue','plum4','lightgoldenrod','lawngreen','forestgreen','dimgray','deeppink','red2','paleturquoise1','palevioletred',
  '#E0CA3C','#989FCE','#261132','#F37748','#AFFC41','#347FC4','#0B5351','#00A9A5',
  '#4E8098','#D81159','#AB3428','#01BAEF','#F42C04','#610345','#044B7F',
  'orchid4','purple4','plum1','olivedrab2','slateblue','mediumvioletred','sienna','orange','seagreen','lightseagreen','mediumpurple4','#F8F4A6',
  'white','yellow','#A88FAC','#D2D68D','#CBD4C2','red','blue','purple','magenta'))

ld <- function(f) { e <- new.env(); load(f, envir = e); get(ls(e)[1], envir = e) }
md   <- ld(file.path(src, "clustering vs features experiment/meta.data.feature.depth.experiment.2025-07-12.Robj"))
emb  <- ld(file.path(src, "clustering vs features experiment/embedding.data.feature.depth.experiment.2025-07-12.Robj"))
rawe <- ld(file.path(src, "clustering vs features ELSA/clustering.by.nfeatures.raw.elsa.values.Robj"))
tre  <- ld(file.path(src, "clustering vs features ELSA/clustering.by.nfeatures.mean.elsa.trends.Robj"))
stopifnot(length(rawe) == length(FEAT_NUM), length(tre) == length(FEAT_NUM))

## variance ranking used to order the sweep (Sup7: vf_vst_counts_variance.standardized)
obj <- ld(file.path(src, "clustering vs features experiment/object.feature.depth.experiment.2025-07-12.Robj"))
vm <- obj@assays$m.est.min@meta.data
vcol <- grep("variance.standardized$", colnames(vm), value = TRUE)[1]
var_sc <- sort(vm[[vcol]], decreasing = TRUE); rm(obj); invisible(gc())

## ELSA rasters: values are in raster cell order (rows from the top, x increasing within rows)
UX <- sort(unique(md$x)); UY <- sort(unique(md$y))
rgrid <- data.table(x = rep(UX, times = length(UY)), y = rep(rev(UY), each = length(UX)))

TXT <- 8
dark <- theme_void(base_size = TXT, base_family = "Liberation Sans") +
  theme(plot.background = element_rect(fill = "black", colour = NA),
        panel.background = element_rect(fill = "black", colour = NA),
        plot.title = element_text(colour = "grey92", size = TXT, hjust = 0.5,
                                  margin = margin(1, 0, 2, 0)),
        plot.margin = margin(2, 2, 2, 2), legend.position = "none")

umap_panel <- function(n) {
  e <- emb[[paste0("UMAP", n, "feature")]]@cell.embeddings
  cl <- md[[paste0(n, "featureClusters")]]
  d <- data.table(u1 = e[, 1], u2 = e[, 2], cl = cl)
  ggplot(d, aes(u1, u2, colour = cl)) + geom_point(size = 0.12, stroke = 0) +
    scale_colour_manual(values = cols.global) + coord_fixed() + dark + ggtitle(lab_depth(n))
}
domain_panel <- function(n) {
  d <- data.table(x = md$x, y = md$y, cl = md[[paste0(n, "featureClusters")]])
  ggplot(d, aes(x, y, fill = cl)) + geom_raster() +
    scale_fill_manual(values = cols.global) + coord_fixed(ratio = 1) + dark +   # the aspect fix
    ggtitle(sprintf("%d domains", length(unique(d$cl)))) +
    ## the sub-region rendered in 3D below (Sup8: x 18000-20000, y 16000-18000), as published
    (if (n == max(DEPTHS)) annotate("rect", xmin = 18000, xmax = 20000, ymin = 16000, ymax = 18000,
                                    fill = NA, colour = "white", linetype = "dashed", linewidth = 0.4))
}
elsa_panel <- function(n, legend = FALSE) {
  v <- rawe[[match(n, FEAT_NUM)]][[1]]                 # d = 50, as in the published figure
  d <- cbind(rgrid, v = v)
  p <- ggplot(d, aes(x, y, fill = v)) + geom_raster() +
    scale_fill_viridis_c(option = "A", limits = c(0, 1), breaks = c(0, 0.5, 1), name = "ELSA") +
    coord_fixed(ratio = 1) + dark
  if (legend) p <- p + theme(legend.position = "right", legend.key.width = unit(0.08, "in"),
                             legend.key.height = unit(0.22, "in"),
                             legend.text = element_text(colour = "grey85", size = TXT - 1),
                             legend.title = element_text(colour = "grey90", size = TXT - 1))
  p
}

## summary curves
ord <- rbindlist(lapply(seq_along(FEAT_NUM), function(i)
  data.table(depth = FEAT_NUM[i], radius = tre[[i]]$d[1:3], elsa = tre[[i]]$elsa[1:3])))
ord_dom <- data.table(depth = FEAT_NUM,
                      domains = sapply(FEAT_NUM, function(n) length(unique(md[[paste0(n, "featureClusters")]]))))
f <- unlist(lapply(c("results-seurat-d10", "results-seurat-d50", "results-seurat", "results-seurat-350",
                     "results-seurat-d650", "results-seurat-d1000"),
                   function(k) list.files(file.path(sw, k), "csv$", full.names = TRUE)))
rnd <- rbindlist(lapply(f, fread), fill = TRUE)[arm == "random"]
rs <- rnd[, .(m = mean(elsa_d50), s = sd(elsa_d50), dm = mean(domains), ds = sd(domains), n = .N), by = depth][order(depth)]
cat("random-subset summary:\n"); print(rs)

light <- theme_classic(base_size = TXT, base_family = "Liberation Sans") +
  theme(plot.background = element_rect(fill = "black", colour = NA),
        panel.background = element_rect(fill = "black", colour = NA),
        axis.text = element_text(colour = "grey85", size = TXT), axis.title = element_text(colour = "grey90", size = TXT),
        axis.line = element_line(colour = "grey70", linewidth = 0.3), axis.ticks = element_line(colour = "grey70", linewidth = 0.3),
        plot.title = element_text(colour = "grey92", size = TXT), legend.background = element_rect(fill = "black"),
        legend.text = element_text(colour = "grey85", size = TXT - 1), legend.title = element_text(colour = "grey90", size = TXT - 1),
        legend.key = element_rect(fill = "black"), plot.margin = margin(3, 6, 3, 3))
vl <- geom_vline(xintercept = DEPTHS, colour = "grey45", linewidth = 0.25, linetype = "dotted")
ORD <- "#E8C547"; RND <- "#4FB3BF"
pA <- ggplot() + vl +
  geom_line(data = ord[radius == 50], aes(depth, elsa, colour = "variance order"), linewidth = 0.5) +
  geom_point(data = ord[radius == 50], aes(depth, elsa, colour = "variance order"), size = 0.6) +
  geom_errorbar(data = rs, aes(depth, ymin = m - s, ymax = m + s), colour = RND, width = 0, linewidth = 0.4) +
  geom_line(data = rs, aes(depth, m, colour = "random subsets"), linewidth = 0.5) +
  geom_point(data = rs, aes(depth, m, colour = "random subsets"), size = 0.8) +
  scale_colour_manual(values = c("variance order" = ORD, "random subsets" = RND), name = NULL) +
  labs(x = "interaction fields entering", y = "mean ELSA, radius 50") + light +
  theme(legend.position = c(0.68, 0.86), legend.key.size = unit(0.12, "in"))
pB <- ggplot(data.table(rank = seq_along(var_sc), v = var_sc), aes(rank, v)) + vl +
  geom_col(fill = "grey75", width = 1) +
  labs(x = "rank of interaction field", y = "standardized variance") + light
pC <- ggplot() + vl +
  geom_line(data = ord_dom, aes(depth, domains), colour = ORD, linewidth = 0.5) +
  geom_point(data = ord_dom, aes(depth, domains), colour = ORD, size = 0.6) +
  geom_errorbar(data = rs, aes(depth, ymin = dm - ds, ymax = dm + ds), colour = RND, width = 0, linewidth = 0.4) +
  geom_line(data = rs, aes(depth, dm), colour = RND, linewidth = 0.5) +
  geom_point(data = rs, aes(depth, dm), colour = RND, size = 0.8) +
  labs(x = "interaction fields entering", y = "number of domains") + light

## 3D surfaces from the published figure
img <- readPNG(pub_img)
strip <- img[4650:dim(img)[1], , 1:3]
strip_g <- wrap_elements(full = grid::rasterGrob(strip, interpolate = TRUE))

rowlab <- function(s) wrap_elements(full = grid::textGrob(s, rot = 90, gp = grid::gpar(col = "grey92", fontsize = TXT, fontfamily = "Liberation Sans")))
blank_lab <- function() wrap_elements(full = grid::rectGrob(gp = grid::gpar(fill = "black", col = NA)))

top <- wrap_plots(c(list(rowlab("UMAP")),   lapply(DEPTHS, umap_panel),
                    list(rowlab("domains")), lapply(DEPTHS, domain_panel),
                    list(rowlab("ELSA")),    lapply(DEPTHS, function(n) elsa_panel(n, legend = n == 1481))),
                  ncol = 6, widths = c(0.12, 1, 1, 1, 1, 1), heights = c(1.05, 0.72, 0.72))
summ <- (pA | pB | pC)
W <- 7.5
## strip aspect fixes its height at full width
strip_h <- W * dim(strip)[1] / dim(strip)[2]
H_top <- 3.75; H_sum <- 1.65
H <- H_top + H_sum + strip_h
cat(sprintf("figure %.2f x %.2f in (PLOS max 7.5 x 8.75)\n", W, H))
stopifnot(H <= 8.75)
fig <- wrap_plots(top, summ, strip_g, ncol = 1, heights = c(H_top, H_sum, strip_h)) &
  theme(plot.background = element_rect(fill = "black", colour = NA))
ggsave(file.path(odir, "Fig6-composite.png"), fig, width = W, height = H, dpi = 600, bg = "black")
cat("wrote Fig6-composite.png\n")
