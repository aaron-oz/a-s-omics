## Supporting figure for the plug-in bias of the common-minimum fields (jensen-bias.R).
## Run: distrobox enter emacs-r -- Rscript jensen-bias-figure.R
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
odir <- "/var/home/aoz/Dropbox/genetics/a-s-omics/code-spatial-smoothing/uncertainty"
fdir <- "/var/home/aoz/Dropbox/genetics/a-s-omics/Submission_PLOSOne/submission-2/figures-regenerated"
B <- readRDS(file.path(odir, "jensen-bias-bins.rds")); P <- fread(file.path(odir, "jensen-bias-pairs.csv"))
TXT <- 8
th <- theme_classic(base_size = TXT, base_family = "Liberation Sans") +
  theme(axis.text = element_text(colour = "grey25", size = TXT), axis.title = element_text(size = TXT),
        plot.title = element_text(size = TXT, face = "bold"), axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3), legend.text = element_text(size = TXT - 1),
        legend.title = element_text(size = TXT - 1), plot.margin = margin(3, 5, 3, 3))
## A: where, in the tissue, the bias is large for many pairs (90th percentile over the 1,477 pairs)
pA <- ggplot(B$bins, aes(x, y, fill = 100 * p90_relbias)) + geom_raster() + coord_fixed() +
  scale_fill_gradient(low = "#F1F4F9", high = "#123B6D", name = "%") +
  labs(title = "A  overstatement, 90th pct of pairs") + theme_void(base_size = TXT, base_family = "Liberation Sans") +
  theme(plot.title = element_text(size = TXT, face = "bold"), legend.text = element_text(size = TXT - 1),
        legend.title = element_text(size = TXT - 1), legend.key.width = unit(0.08, "in"),
        legend.key.height = unit(0.2, "in"), plot.margin = margin(3, 3, 3, 3))
## B: the bias concentrates where neither partner is clearly limiting
S <- B$scatter[zgap <= 4]
pB <- ggplot(S, aes(zgap, 100 * rb_med)) + geom_bin_2d(binwidth = c(0.1, 5)) +
  scale_fill_gradient(low = "#DCE3EE", high = "#123B6D", trans = "log10", name = "bins") +
    labs(title = "B  bias vs partner separation", x = "standardized gap between partners",
       y = "overstatement of median (%)") + th + theme(legend.key.width = unit(0.08, "in"))
## C: spatial pattern preserved? per-pair correlation, published vs corrected field
H <- rbind(data.table(r = P$cor_pattern_med,  target = "vs posterior median"),
           data.table(r = P$cor_pattern_mean, target = "vs posterior mean"))
pC <- ggplot(H, aes(r)) + geom_histogram(binwidth = 0.005, fill = "#4A6A93", colour = NA, boundary = 1) +
  facet_wrap(~target, ncol = 1, scales = "free_y") + coord_cartesian(xlim = c(0.85, 1)) +
  labs(title = "C  pattern, per pair", x = "Pearson correlation over bins", y = "pairs") +
  th + theme(strip.background = element_blank(), strip.text = element_text(size = TXT, hjust = 0))
fig <- pA + pB + pC + plot_layout(widths = c(1.25, 1, 0.9))
ggsave(file.path(fdir, "S11_Fig_plugin-bias.png"), fig, width = 7.5, height = 2.6, dpi = 600, bg = "white")
cat("median of per-bin p90 bias:", round(median(B$bins$p90_relbias), 3), " max:", round(max(B$bins$p90_relbias), 3), "\n")
