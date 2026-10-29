## Main-text figures for the uncertainty results, at PLOS size: 7.5 in wide, Liberation Sans,
## no text below 8 pt, 400 dpi. Reads the pilot outputs (the showcase fields, fixed before their
## posteriors were computed). Main text shows the three randomly drawn pairs and Cd44; every
## showcase pair is in the SI version. Files are named by content, not figure number.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript make-main-figs.R
## 2026-10-29 (AOZ review): region figure gains a map of the frozen regions and transects and a
## common ratio axis across its panels; transect figure gains an inset locating the transect;
## whole-field figure gains its color key and a color bar for the plug-in log ratio.
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(ragg)})
pd <- "../pilot"; FONT <- "Liberation Sans"; PT <- 8
th <- theme_bw(base_size = PT, base_family = FONT) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", strip.background = element_rect(fill = "grey92"),
        plot.title = element_text(size = PT, face = "bold"), legend.text = element_text(size = PT), axis.text = element_text(size = PT),
        strip.text = element_text(size = PT))
cat4 <- c("looks different, credible" = "#0072B2", "looks different, not credible" = "#D55E00",
          "looks alike, credible" = "#009E73", "looks alike, not credible" = "#7A7A7A")
shp4 <- setNames(c(16, 17, 15, 1), names(cat4))
lab <- function(ld, cr) factor(ifelse(ld, ifelse(cr, names(cat4)[1], names(cat4)[2]), ifelse(cr, names(cat4)[3], names(cat4)[4])), levels = names(cat4))
nice <- c("box-heart-1" = "heart box 1", "box-heart-2" = "heart box 2", "box-liver" = "liver box", "box-cord-1" = "cord box 1",
          "box-cord-2" = "cord box 2", "anat-heart" = "heart", "anat-liver" = "liver", "anat-lung" = "lung",
          "anat-spinal-cord" = "spinal cord", "anat-thymus" = "thymus")
OD <- "/var/home/aoz/Dropbox/genetics/a-s-omics/Submission_PLOSOne/submission-2/figures-regenerated"
save <- function(g, f, w, h) { agg_png(file.path(OD, f), width = w, height = h, units = "in", res = 400); print(g); invisible(dev.off()) }

## ---- map of the frozen regions and transects (panel A of the region figure; inset of the transect figure)
REPO <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
ann <- fread(file.path(REPO, "data-outputs/mouse-embryo-raw/controls/e16-fov-annotation.csv"))
ann[, `:=`(i = (x - 15025) / 50 + 1, j = (y - 11525) / 50 + 1)]
DS <- readRDS(file.path(REPO, "code-spatial-smoothing/uncertainty/design/design-regions.rds"))
anat5 <- c("Heart", "Liver", "Lung", "Spinal cord", "Thymus")
ann[, ann5 := ifelse(annotation %in% anat5, annotation, "other tissue")]
anatcol <- c(Heart = "#E69F00", Liver = "#56B4E9", Lung = "#009E73", "Spinal cord" = "#CC79A7", Thymus = "#F0E442", "other tissue" = "grey85")
boxes <- copy(DS$boxes)[, lab := c("heart box 1", "heart box 2", "liver box", "cord box 1", "cord box 2")]
slname <- c("slice-cord-to-muscle" = "cord to muscle", "slice-heart-to-liver" = "heart to liver",
            "slice-along-cord" = "along the cord", "slice-thymus-to-heart" = "thymus to heart")
sl <- rbindlist(lapply(names(DS$slice_bins), function(n) { d <- as.data.table(DS$slice_bins[[n]])[order(pos)]; data.table(slice = slname[n], i = d$i, j = d$j) }))
region_map <- function(which_slices = names(slname), show_boxes = TRUE, base = PT, legend = TRUE) {
  g <- ggplot() + geom_raster(data = ann, aes(i, j, fill = ann5)) + coord_fixed(expand = FALSE) +
    scale_fill_manual(values = anatcol, name = NULL, breaks = anat5) +
    theme_void(base_size = base, base_family = FONT) +
    theme(legend.position = if (legend) "right" else "none", legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = base), plot.title = element_text(size = base, face = "bold", hjust = 0))
  d <- sl[slice %in% slname[which_slices]]
  g <- g + geom_path(data = d, aes(i, j, group = slice), color = "black", linewidth = 0.5)
  if (show_boxes) g <- g + geom_rect(data = boxes, aes(xmin = i0 - 0.5, xmax = i1 + 0.5, ymin = j0 - 0.5, ymax = j1 + 0.5),
                                      fill = NA, color = "black", linewidth = 0.5)
  g
}
## labels for the full map (placed by hand, in bin units, clear of the features they name)
box_lab <- boxes[, .(lab, i = (i0 + i1) / 2, j = j1 + 5)]
sl_lab <- data.table(lab = c("cord to muscle", "heart to liver", "along the cord", "thymus to heart"),
                     i = c(57, 199, 153, 42), j = c(40, 90, 6, 79), hjust = c(0, 1, 0, 0))
map_full <- region_map() + ggtitle("A  Frozen regions (boxes and annotated tissues) and transects (lines)") +
  geom_text(data = box_lab, aes(i, j, label = lab), size = PT / .pt, family = FONT) +
  geom_text(data = sl_lab, aes(i, j, label = lab, hjust = hjust), size = PT / .pt, family = FONT)

## ---- region comparisons
regfig <- function(fields, file, ncol, h, map = TRUE, title_b = "B  Ratio of regional averages") {
  r <- fread(file.path(pd, "results-shuffle-regional.csv"))[field %in% fields]
  r[, field := factor(field, levels = fields)]
  r[, cmp := factor(paste(nice[r1], "/", nice[r2]), levels = rev(unique(paste(nice[r1], "/", nice[r2]))))]
  r[, cat := lab(looks_diff, credible)][, y := as.numeric(cmp)]
  g <- ggplot(r) + geom_vline(xintercept = 1, color = "grey50", linewidth = 0.3) +
    annotate("rect", xmin = 1 / 1.25, xmax = 1.25, ymin = -Inf, ymax = Inf, fill = "grey88") +
    geom_segment(aes(x = lo, xend = hi, y = y + 0.12, yend = y + 0.12, color = cat), linewidth = 0.6) +
    geom_segment(aes(x = dec_lo, xend = dec_hi, y = y - 0.22, yend = y - 0.22), color = "black", linewidth = 0.35) +
    geom_point(aes(x = post_ratio, y = y + 0.12, color = cat, shape = cat), size = 1.6) +
    geom_point(aes(x = plug_ratio, y = y + 0.12), shape = 4, size = 1.6) +
    scale_y_continuous(breaks = seq_along(levels(r$cmp)), labels = levels(r$cmp), expand = expansion(add = 0.6)) +
    ## one ratio axis for every panel, so interval widths are comparable across fields
    scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 5), labels = c("0.2", "0.5", "1", "2", "5"), limits = c(0.15, 8)) +
    scale_color_manual(values = cat4, name = NULL, drop = FALSE) + scale_shape_manual(values = shp4, name = NULL, drop = FALSE) +
    facet_wrap(~ field, ncol = ncol) + th + guides(color = guide_legend(nrow = 2)) +
    labs(title = if (map) title_b else NULL, x = "ratio of regional averages, first region / second region (log scale)", y = NULL)
  if (map) g <- (map_full / g) + plot_layout(heights = c(1.45, 3.2))
  save(g, file, 7.5, h)
}
regfig(c("CALR-ITGA2B", "NID1-ITGB1", "ADAM12-ITGB1", "Cd44"), "fig-regions-main.png", 2, 8.6)
regfig(c("CALR-ITGA2B", "NID1-ITGB1", "ADAM12-ITGB1", "WNT4-FZD6", "RSPO1-LGR6", "DHH-PTCH1", "VEGFA-KDR", "Cd44"),
       "fig-regions-all-SI.png", 2, 11, map = FALSE)

## ---- transects: the drawn transect pair and Cd44, both on the heart-to-liver transect
s <- fread(file.path(pd, "results-slices.csv"))[slice == "slice-heart-to-liver" & field %in% c("ADAM12-ITGB1", "Cd44")]
s[, field := factor(field, levels = c("ADAM12-ITGB1", "Cd44"))]
s[, tissue := ifelse(annotation %in% c("Heart", "Liver"), annotation, "other")]
tr <- function(f) {
  d <- s[field == f]; lo <- min(d$sim_lo); hi <- max(d$sim_hi); strip <- lo - 0.06 * (hi - lo)
  ggplot(d, aes(pos)) + geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    geom_ribbon(aes(ymin = sim_lo, ymax = sim_hi), fill = "#0072B2", alpha = 0.18) +
    geom_ribbon(aes(ymin = pw_lo, ymax = pw_hi), fill = "#0072B2", alpha = 0.30) +
    geom_line(aes(y = med), color = "#0072B2", linewidth = 0.5) +
    geom_line(aes(y = plug_dev), color = "black", linewidth = 0.35, linetype = "22") +
    geom_tile(aes(y = strip, fill = tissue), height = 0.05 * (hi - lo)) +
    geom_point(data = d[exc_pos | exc_neg], aes(y = strip - 0.06 * (hi - lo)), shape = 124, size = 1.8, color = "#D55E00") +
    scale_fill_manual(values = c(Heart = "#E69F00", Liver = "#56B4E9", other = "grey80"), name = "tissue") +
    th + labs(title = f, x = "position along the heart-to-liver transect (bins)", y = "log(field / transect average)")
}
inset <- region_map("slice-heart-to-liver", show_boxes = FALSE, legend = FALSE) + guides(fill = "none") +
  theme(plot.background = element_rect(fill = "white", color = "grey40", linewidth = 0.3))
tr_a <- tr("ADAM12-ITGB1") + inset_element(inset, left = 0.70, bottom = 0.74, right = 0.995, top = 0.995, align_to = "panel")
save((tr_a | tr("Cd44")) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
     "fig-transects-main.png", 7.5, 3.4)

## ---- whole fields
W <- readRDS("/var/home/aoz/data/a-s-omics-uncertainty/wholefield-additions.rds")
xy <- W$xy[, .(i = (x - 15025) / 50 + 1, j = (y - 11525) / 50 + 1)]
cols <- c(above = "#D55E00", below = "#0072B2", neither = "grey92")
lab3 <- function(p, n) factor(ifelse(p, "above", ifelse(n, "below", "neither")), c("above", "below", "neither"))
mth <- theme_void(base_size = PT, base_family = FONT) + theme(plot.title = element_text(size = PT, hjust = 0), legend.position = "none")
row_for <- function(nm, short) {
  b <- W$maps[[paste(nm, "bin")]]; k <- W$maps[[paste(nm, "blk")]]
  d <- cbind(xy, plug = pmax(pmin(b$plug, 2), -2), pw = lab3(b$J$pw_pos, b$J$pw_neg), ex = lab3(b$J$pos, b$J$neg),
             exb = lab3(k$J$pos[W$blk], k$J$neg[W$blk]))
  m <- function(v, t) ggplot(d, aes(i, j, fill = .data[[v]])) + geom_raster() + coord_fixed(expand = FALSE) + mth +
    scale_fill_manual(values = cols, drop = FALSE) + labs(title = t)
  pa <- ggplot(d, aes(i, j, fill = plug)) + geom_raster() + coord_fixed(expand = FALSE) + mth +
    scale_fill_gradient2(low = "#0072B2", mid = "grey92", high = "#D55E00", midpoint = 0, limits = c(-2, 2)) +
    labs(title = paste0(short, "\nplug-in log ratio"))
  pa | m("pw", sprintf("\npointwise: %s bins", format(sum(b$J$pw_pos | b$J$pw_neg), big.mark = ","))) |
    m("ex", sprintf("\njoint, bins: %s", format(sum(b$J$pos | b$J$neg), big.mark = ","))) |
    m("exb", sprintf("\njoint, 5 x 5 blocks: %d of 1,040", sum(k$J$pos | k$J$neg)))
}
## keys drawn as small plots (a shared legend cannot be collected across the two fill scales)
cbar <- ggplot(data.frame(v = seq(-2, 2, length.out = 201)), aes(v, 1, fill = v)) + geom_raster() +
  scale_fill_gradient2(low = "#0072B2", mid = "grey92", high = "#D55E00", midpoint = 0, limits = c(-2, 2), guide = "none") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), labels = c("-2 and below", "-1", "0", "1", "2 and above"), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.5, 1.5)) + theme_void(base_size = PT, base_family = FONT) +
  theme(axis.text.x = element_text(size = PT, margin = margin(t = 2)), plot.title = element_text(size = PT, hjust = 0.5)) +
  ggtitle("plug-in log ratio")
key <- ggplot(data.frame(x = c(1, 2.6), f = factor(c("above", "below"), names(cols)), lab = c("credibly above", "credibly below"))) +
  geom_tile(aes(x, 1, fill = f), width = 0.4, height = 0.6) + geom_text(aes(x + 0.25, 1, label = lab), hjust = 0, size = PT / .pt, family = FONT) +
  scale_fill_manual(values = cols, guide = "none") + coord_cartesian(xlim = c(0.7, 4.2), ylim = c(0.5, 1.5)) +
  theme_void(base_size = PT, base_family = FONT) + theme(plot.title = element_text(size = PT, hjust = 0.5)) +
  ggtitle("posterior calls (grey: neither)")
g <- row_for("CALM1-KCNQ3 vs CALM1-PDE1A", "A  CALM1-KCNQ3 / CALM1-PDE1A") / row_for("Cd44 vs its tissue average", "B  Cd44 / its average") /
  row_for("PKM-CD44 vs its tissue average", "C  PKM-CD44 / its average") /
  ((plot_spacer() | cbar | plot_spacer() | key | plot_spacer()) + plot_layout(widths = c(0.3, 1.6, 0.3, 1.6, 0.3)))
g <- g + plot_layout(heights = c(1, 1, 1, 0.33))
save(g, "fig-wholefield-main.png", 7.5, 5.1)
cat("done\n")

## ---- SI: cross-field regional comparison (added to S12 Fig under the all-pairs region figure)
x <- fread(file.path(pd, "results-crossfield-regional.csv"))
x <- merge(x, fread(file.path(pd, "results-shuffle-crossfield.csv"))[, .(region, dec_lo, dec_hi)], by = "region")
x[, cat := lab(looks_diff, credible)][, reg := factor(nice[region], levels = rev(nice))]
gx <- ggplot(x, aes(y = reg)) + geom_vline(xintercept = 1, color = "grey50", linewidth = 0.3) +
  annotate("rect", xmin = 1 / 1.25, xmax = 1.25, ymin = -Inf, ymax = Inf, fill = "grey88") +
  geom_segment(aes(x = lo, xend = hi, yend = reg, color = cat), linewidth = 0.6) +
  geom_point(aes(x = post_ratio, color = cat, shape = cat), size = 1.6) + geom_point(aes(x = plug_ratio), shape = 4, size = 1.6) +
  scale_x_log10() + scale_color_manual(values = cat4, name = NULL, drop = FALSE) + scale_shape_manual(values = shp4, name = NULL, drop = FALSE) +
  th + guides(color = guide_legend(nrow = 2)) +
  labs(title = "CALM1-KCNQ3 / CALM1-PDE1A, each divided by its own average, within each region", x = "ratio (log scale)", y = NULL)
save(gx, "fig-crossfield-SI.png", 7.5, 3.2)
