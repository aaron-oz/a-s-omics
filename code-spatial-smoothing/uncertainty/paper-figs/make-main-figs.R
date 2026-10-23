## Main-text figures for the uncertainty results, at PLOS size: 7.5 in wide, Liberation Sans,
## no text below 8 pt, 400 dpi. Reads the pilot outputs (the showcase fields, fixed before their
## posteriors were computed). Main text shows the three randomly drawn pairs and Cd44; every
## showcase pair is in the SI version. Files are named by content, not figure number.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript make-main-figs.R
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

## ---- region comparisons
regfig <- function(fields, file, ncol, h) {
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
    scale_x_log10() + scale_color_manual(values = cat4, name = NULL, drop = FALSE) + scale_shape_manual(values = shp4, name = NULL, drop = FALSE) +
    facet_wrap(~ field, ncol = ncol, scales = "free_x") + th + guides(color = guide_legend(nrow = 2)) +
    labs(x = "ratio of regional averages, first region / second region (log scale)", y = NULL)
  save(g, file, 7.5, h)
}
regfig(c("CALR-ITGA2B", "NID1-ITGB1", "ADAM12-ITGB1", "Cd44"), "fig-regions-main.png", 2, 6.2)
regfig(c("CALR-ITGA2B", "NID1-ITGB1", "ADAM12-ITGB1", "WNT4-FZD6", "RSPO1-LGR6", "DHH-PTCH1", "VEGFA-KDR", "Cd44"),
       "fig-regions-all-SI.png", 2, 11)

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
save((tr("ADAM12-ITGB1") | tr("Cd44")) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
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
key <- ggplot(data.frame(x = 1:3, f = factor(names(cols), names(cols)))) + geom_tile(aes(x, 1, fill = f)) +
  scale_fill_manual(values = cols, name = "credibly") + theme_void(base_size = PT, base_family = FONT) + theme(legend.position = "bottom")
g <- row_for("CALM1-KCNQ3 vs CALM1-PDE1A", "A  CALM1-KCNQ3 / CALM1-PDE1A") / row_for("Cd44 vs its tissue average", "B  Cd44 / its average") /
  row_for("PKM-CD44 vs its tissue average", "C  PKM-CD44 / its average")
save(g, "fig-wholefield-main.png", 7.5, 4.6)
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
