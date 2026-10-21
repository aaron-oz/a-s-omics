## Draft pilot figures (Phase 3, step 6). Reads the pilot-run.R outputs.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript pilot-figures.R
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
cat4 <- c("looks different, credible" = "#0072B2", "looks different, NOT credible" = "#D55E00",
          "looks alike, credible" = "#009E73", "looks alike, not credible" = "#7A7A7A")
shp4 <- c(16, 17, 15, 1); names(shp4) <- names(cat4)
lab <- function(ld, cr) factor(ifelse(ld, ifelse(cr, names(cat4)[1], names(cat4)[2]),
                                      ifelse(cr, names(cat4)[3], names(cat4)[4])), levels = names(cat4))
th <- theme_bw(9) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")

## ---- Fig P1: regional, same field
r <- fread("results-shuffle-regional.csv")   # results-regional.csv plus the decorrelated intervals
r[, cat := lab(looks_diff, credible)][, cmp := paste(r1, "vs", r2)]
r[, wil := ifelse(wilcox_p < 0.05, "Wilcoxon p < 0.05", "Wilcoxon p >= 0.05")]
g1 <- ggplot(r, aes(y = cmp)) +
  geom_vline(xintercept = 1, color = "grey50", linewidth = 0.3) +
  annotate("rect", xmin = 1 / 1.25, xmax = 1.25, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
  geom_segment(aes(x = lo, xend = hi, yend = cmp, color = cat), linewidth = 0.6) +
  geom_segment(aes(x = dec_lo, xend = dec_hi, y = as.numeric(factor(cmp)) - 0.3, yend = as.numeric(factor(cmp)) - 0.3),
               color = "black", linewidth = 0.4) +
  geom_point(aes(x = post_ratio, color = cat, shape = cat), size = 1.8) +
  geom_point(aes(x = plug_ratio), shape = 4, size = 1.6) +
  scale_x_log10() + scale_color_manual(values = cat4, name = NULL) + scale_shape_manual(values = shp4, name = NULL) +
  facet_wrap(~ field, ncol = 4) + th + labs(x = "ratio of regional averages, region 1 / region 2 (log scale)", y = NULL,
    title = "Regional comparisons, same field: posterior 95% interval (colored) vs plug-in ratio (x)",
    subtitle = "Gray band: plug-in ratio within 1/1.25 to 1.25 ('looks alike'). Credible = 95% interval excludes 1.\nThin black bar below each row: 95% interval if spatial correlation is ignored (draws shuffled across bins).")
ggsave("fig-P1-regional.png", g1, width = 13, height = 9, dpi = 120)

## ---- Fig P2: regional, cross field
x <- fread("results-crossfield-regional.csv")[, cat := lab(looks_diff, credible)]
g2 <- ggplot(x, aes(y = region)) + geom_vline(xintercept = 1, color = "grey50", linewidth = 0.3) +
  annotate("rect", xmin = 1 / 1.25, xmax = 1.25, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
  geom_segment(aes(x = lo, xend = hi, yend = region, color = cat), linewidth = 0.8) +
  geom_point(aes(x = post_ratio, color = cat, shape = cat), size = 2) +
  geom_point(aes(x = plug_ratio), shape = 4, size = 2) +
  geom_text(aes(x = 60, label = sprintf("Wilcoxon signed-rank p = %.1g", wilcox_signed_p)), size = 2.6, hjust = 1) +
  scale_x_log10(limits = c(0.05, 60)) + scale_color_manual(values = cat4, name = NULL) + scale_shape_manual(values = shp4, name = NULL) +
  th + labs(x = "shape-normalized ratio, CALM1-KCNQ3 / CALM1-PDE1A (log scale)", y = NULL,
            title = "Cross-field comparison within each region (each field divided by its own tissue average)")
ggsave("fig-P2-crossfield.png", g2, width = 9, height = 4.5, dpi = 120)

## ---- Fig P3: slices (the frozen slice slot, ADAM12-ITGB1, plus Cd44 and NID1-ITGB1 on two slices)
s <- fread("results-slices.csv")
sel <- s[(field == "ADAM12-ITGB1" & slice == "slice-heart-to-liver") | (field == "Cd44" & slice == "slice-heart-to-liver") |
         (field == "NID1-ITGB1" & slice == "slice-cord-to-muscle") | (field == "Cd44" & slice == "slice-along-cord")]
sel[, panel := paste(field, slice, sep = ", ")]
g3 <- ggplot(sel, aes(pos)) + geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_ribbon(aes(ymin = sim_lo, ymax = sim_hi), fill = "#0072B2", alpha = 0.18) +
  geom_ribbon(aes(ymin = pw_lo, ymax = pw_hi), fill = "#0072B2", alpha = 0.30) +
  geom_line(aes(y = med), color = "#0072B2", linewidth = 0.6) +
  geom_line(aes(y = plug_dev), color = "black", linewidth = 0.4, linetype = "22") +
  geom_rug(data = sel[exc_pos | exc_neg], sides = "b", color = "#D55E00", linewidth = 0.6) +
  facet_wrap(~ panel, scales = "free", ncol = 2) + th +
  labs(x = "position along slice (bins)", y = "log(field / slice average)",
       title = "Slices: plug-in (dashed) vs posterior median (line), pointwise 95% (dark band), simultaneous 95% (light band)",
       subtitle = "Orange ticks: positions in the joint 95% excursion set (credibly above or below the slice average)")
ggsave("fig-P3-slices.png", g3, width = 11, height = 7, dpi = 120)

## ---- Fig P4: whole fields, three comparisons (rows) x four panels (columns), from additions-run.R
A <- readRDS("/var/home/aoz/data/a-s-omics-uncertainty/wholefield-additions.rds")
xy <- A$xy[, .(i = (x - 15025) / 50 + 1, j = (y - 11525) / 50 + 1)]
cols <- c(above = "#D55E00", below = "#0072B2", neither = "grey92")
lab3 <- function(p, n) factor(ifelse(p, "above", ifelse(n, "below", "neither")), c("above", "below", "neither"))
base <- list(coord_fixed(expand = FALSE), theme_void(8), theme(legend.position = "none", plot.title = element_text(size = 7.5)))
row_for <- function(nm) {
  b <- A$maps[[paste(nm, "bin")]]; k <- A$maps[[paste(nm, "blk")]]
  kb <- k$J; kb <- lapply(kb[c("pos", "neg")], function(v) v[A$blk])   # block results painted onto bins
  d <- cbind(xy, plug = pmax(pmin(b$plug, 3), -3), pw = lab3(b$J$pw_pos, b$J$pw_neg),
             ex = lab3(b$J$pos, b$J$neg), exb = lab3(kb$pos, kb$neg))
  pa <- ggplot(d, aes(i, j, fill = plug)) + geom_raster() + base + theme(legend.position = "right") +
    scale_fill_gradient2(low = "#0072B2", mid = "grey92", high = "#D55E00", midpoint = 0, name = "log ratio") +
    labs(title = paste0(nm, "\nplug-in (what the map shows)"))
  f <- function(v, t) ggplot(d, aes(i, j, fill = .data[[v]])) + geom_raster() + base + scale_fill_manual(values = cols) + labs(title = t)
  pa | f("pw", sprintf("\npointwise P >= 0.95: %d bins", sum(b$J$pw_pos | b$J$pw_neg))) |
    f("ex", sprintf("\njoint 95%%, bins: %d of 26,000", sum(b$J$pos | b$J$neg))) |
    f("exb", sprintf("\njoint 95%%, 5x5 blocks: %d of 1,040", sum(k$J$pos | k$J$neg)))
}
g4 <- row_for("CALM1-KCNQ3 vs CALM1-PDE1A") / row_for("Cd44 vs its tissue average") / row_for("PKM-CD44 vs its tissue average") +
  plot_annotation(title = "Whole fields: orange = credibly above, blue = credibly below (one field vs another, or a field vs its own tissue average)")
ggsave("fig-P4-wholefield.png", g4, width = 15, height = 9.5, dpi = 120)
