## Graphical abstract for the uncertainty reframing (2026-10-23). Uses ADAM12-ITGB1, a pair drawn
## at random under the frozen rule, compared between the two heart boxes; and the systematic
## counts (Table 2). Replaces graphical-abstract.R (domains), kept for the record.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript graphical-abstract-uq.R
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(ragg)})
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"; FONT <- "Liberation Sans"
mn <- readRDS(file.path(repo, "code-spatial-smoothing/controls/sweep/m-est-min-matrix.rds"))
keep <- !grepl("conflicted", mn$features); M <- mn$counts[keep, ]; rownames(M) <- mn$features[keep]
tc <- readRDS(file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
D <- readRDS(file.path(repo, "code-spatial-smoothing/uncertainty/design/design-regions.rds"))
r <- fread(file.path(repo, "code-spatial-smoothing/uncertainty/pilot/results-shuffle-regional.csv"))[field == "ADAM12-ITGB1" & r1 == "box-heart-1" & r2 == "box-heart-2"]
d <- data.table(i = (tc$x - 15025) / 50 + 1, j = (tc$y - 11525) / 50 + 1, v = M["ADAM12-ITGB1", ] / tc$total.count)
d[, v := pmin(v, quantile(v, 0.99))]
bx <- D$boxes[name %in% c("box-heart-1", "box-heart-2")]
th <- theme_void(base_size = 11, base_family = FONT) + theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
  plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey25"), legend.position = "none", plot.margin = margin(6, 8, 6, 8))
pA <- ggplot(d, aes(i, j, fill = v)) + geom_raster() + scale_fill_viridis_c() + coord_fixed(expand = FALSE) +
  geom_rect(data = bx, aes(xmin = i0 - .5, xmax = i1 + .5, ymin = j0 - .5, ymax = j1 + .5), inherit.aes = FALSE, fill = NA,
            colour = c("white", "#FF6B6B"), linewidth = 0.9) +
  annotate("text", x = (bx$i0 + bx$i1) / 2, y = bx$j1 + 6, label = c("A", "B"), colour = c("white", "#FF6B6B"), family = FONT, fontface = "bold", size = 5) +
  th + labs(title = "The map", subtitle = sprintf("ligand-receptor field is %.1fx higher in A than in B", r$plug_ratio))
iv <- data.table(reading = factor(c("Map read by eye", "Wilcoxon, bins independent", "Posterior, correlation removed", "Joint posterior"),
                                   levels = rev(c("Map read by eye", "Wilcoxon, bins independent", "Posterior, correlation removed", "Joint posterior"))),
                 lo = c(NA, NA, r$dec_lo, r$lo), hi = c(NA, NA, r$dec_hi, r$hi), mid = c(r$plug_ratio, NA, NA, r$post_ratio),
                 verdict = c("different", sprintf("different (p = %.4f)", r$wilcox_p), "different", "not credibly different"))
pB <- ggplot(iv, aes(y = reading)) + geom_vline(xintercept = 1, colour = "grey50") +
  geom_segment(aes(x = lo, xend = hi, yend = reading, colour = reading == "Joint posterior"), linewidth = 2, na.rm = TRUE) +
  geom_point(aes(x = mid), shape = 4, size = 3, stroke = 1.2, na.rm = TRUE) +
  geom_text(aes(x = 3.6, label = verdict, colour = reading == "Joint posterior"), hjust = 0, family = FONT, size = 3.9) +
  scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2")) +
  scale_x_log10(limits = c(0.6, 30), breaks = c(1, 1.5, 2.5)) + theme_minimal(base_size = 11, base_family = FONT) +
  theme(legend.position = "none", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey25"), axis.text.y = element_text(size = 11, colour = "black")) +
  labs(title = "Is A really higher than B?", subtitle = "ratio A / B, with 95% intervals where the reading has one", x = "ratio (log scale)", y = NULL)
sm <- data.table(reading = factor(c("Wilcoxon, bins independent", "Posterior, correlation removed", "Joint posterior"),
                                  levels = rev(c("Wilcoxon, bins independent", "Posterior, correlation removed", "Joint posterior"))),
                 pct = c(90, 73, 46))
pC <- ggplot(sm, aes(pct, reading, fill = reading == "Joint posterior")) + geom_col(width = 0.6) +
  geom_text(aes(label = paste0(pct, "%")), hjust = -0.15, family = FONT, size = 4) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2")) + scale_x_continuous(limits = c(0, 105), breaks = NULL) +
  theme_minimal(base_size = 11, base_family = FONT) + theme(legend.position = "none", panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5), plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey25"),
    axis.text.y = element_text(size = 11, colour = "black")) +
  labs(title = "Across 1,477 ligand-receptor pairs", subtitle = "share of 19,201 region comparisons called different", x = NULL, y = NULL)
g <- (pA | (pB / pC)) + plot_layout(widths = c(1.15, 1)) +
  plot_annotation(title = "Do visible differences in ligand-receptor fields survive their uncertainty?",
                  caption = "One E16.5 mouse embryo section (Stereo-seq); the 19,201 comparisons share regions and pairs, so the shares are not rates. Fields estimate the potential for signaling, not signaling itself.",
                  theme = theme(plot.title = element_text(family = FONT, face = "bold", size = 16, hjust = 0.5),
                                plot.caption = element_text(family = FONT, size = 10, colour = "grey35", hjust = 0.5)))
agg_png("graphical-abstract-uq.png", width = 12, height = 7, units = "in", res = 250); print(g); invisible(dev.off())
cat("ok\n")
