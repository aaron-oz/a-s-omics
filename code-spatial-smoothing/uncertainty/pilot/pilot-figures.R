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
r <- fread("results-regional.csv")
r[, cat := lab(looks_diff, credible)][, cmp := paste(r1, "vs", r2)]
r[, wil := ifelse(wilcox_p < 0.05, "Wilcoxon p < 0.05", "Wilcoxon p >= 0.05")]
g1 <- ggplot(r, aes(y = cmp)) +
  geom_vline(xintercept = 1, color = "grey50", linewidth = 0.3) +
  annotate("rect", xmin = 1 / 1.25, xmax = 1.25, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
  geom_segment(aes(x = lo, xend = hi, yend = cmp, color = cat), linewidth = 0.6) +
  geom_point(aes(x = post_ratio, color = cat, shape = cat), size = 1.8) +
  geom_point(aes(x = plug_ratio), shape = 4, size = 1.6) +
  scale_x_log10() + scale_color_manual(values = cat4, name = NULL) + scale_shape_manual(values = shp4, name = NULL) +
  facet_wrap(~ field, ncol = 4) + th + labs(x = "ratio of regional averages, region 1 / region 2 (log scale)", y = NULL,
    title = "Regional comparisons, same field: posterior 95% interval (colored) vs plug-in ratio (x)",
    subtitle = "Gray band: plug-in ratio within 1/1.25 to 1.25 ('looks alike'). Credible = 95% interval excludes 1.")
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

## ---- Fig P4/P5: whole fields
W <- readRDS("/var/home/aoz/data/a-s-omics-uncertainty/wholefield.rds")
xy <- W$xy[, .(i = (x - 15025) / 50 + 1, j = (y - 11525) / 50 + 1)]
wf_panels <- function(J, plugv, ttl) {
  d <- cbind(xy, plug = pmax(pmin(plugv, 3), -3),
             pw = factor(ifelse(J$pw_pos, "above", ifelse(J$pw_neg, "below", "neither")), c("above", "below", "neither")),
             ex = factor(ifelse(J$pos, "above", ifelse(J$neg, "below", "neither")), c("above", "below", "neither")))
  cols <- c(above = "#D55E00", below = "#0072B2", neither = "grey92")
  base <- list(coord_fixed(expand = FALSE), theme_void(9), theme(legend.position = "bottom"))
  a <- ggplot(d, aes(i, j, fill = plug)) + geom_raster() + base +
    scale_fill_gradient2(low = "#0072B2", mid = "grey92", high = "#D55E00", midpoint = 0, name = "plug-in log ratio") +
    labs(title = "plug-in (what the map shows)")
  b <- ggplot(d, aes(i, j, fill = pw)) + geom_raster() + base + scale_fill_manual(values = cols, name = NULL) +
    labs(title = sprintf("pointwise P >= 0.95: %d above, %d below", sum(J$pw_pos), sum(J$pw_neg)))
  c <- ggplot(d, aes(i, j, fill = ex)) + geom_raster() + base + scale_fill_manual(values = cols, name = NULL) +
    labs(title = sprintf("joint 95%% excursion sets: %d above, %d below", sum(J$pos), sum(J$neg)))
  (a | b | c) + plot_annotation(title = ttl)
}
ggsave("fig-P4-wholefield-cross.png", wf_panels(W$Jx, W$plug_x, "Whole field: log(CALM1-KCNQ3 / CALM1-PDE1A), each shape-normalized"),
       width = 14, height = 4.2, dpi = 120)
ggsave("fig-P5-wholefield-cd44.png", wf_panels(W$Jc, W$plug_c, "Whole field: log(Cd44 / its tissue average)"),
       width = 14, height = 4.2, dpi = 120)
