## Fig 1, project schematic, redrawn in code (2026-10-23) from the original hand-made image
## (plos-port/img/Project Schematic.png) with the same layout, plus the uncertainty step:
## a posterior-density glyph on each estimated quantity and a "Credible differences" box fed by
## the joint posterior. 7.5 in wide; no text below 8 pt. Liberation Sans.
## Run: distrobox enter emacs-r -- Rscript fig1-schematic.R
suppressMessages({library(ggplot2); library(ragg)})
FONT <- "Liberation Sans"; PUR <- "#8E00E6"; RED <- "#E8262B"; ORA <- "#F2921D"; TEA <- "#0B7E98"
pt <- function(p) p / .pt                       # ggplot text size from points
blob <- function(cx, cy, rx, ry, col, n = 45, amax = 0.05) {   # soft glow: stacked translucent ellipses
  do.call(rbind, lapply(seq_len(n), function(k) { f <- k / n; t <- seq(0, 2 * pi, length.out = 120)
    data.frame(x = cx + rx * f * cos(t), y = cy + ry * f * sin(t), g = paste(cx, cy, k), col = col, a = amax) }))
}
B <- rbind(blob(6.9, 4.6, 3.1, 1.5, ORA, amax = 0.016), blob(5.0, 2.0, 4.6, 1.9, TEA, amax = 0.018))
boxes <- data.frame(
  lab = c("Cell states\nin tissues", "L/R mRNA\ndensity", "L/R protein\ndensity", "L-R binding\ndensity", "Cell-to-cell\ncommunication"),
  x = c(5.0, 8.3, 6.8, 3.2, 1.7), y = c(6.35, 4.3, 1.75, 1.75, 4.3),
  est = c(FALSE, TRUE, TRUE, TRUE, FALSE), hd = c(FALSE, TRUE, TRUE, TRUE, TRUE))
W <- 1.75; H <- 0.95
stack <- do.call(rbind, lapply(which(boxes$hd), function(i) do.call(rbind, lapply(4:1, function(k)
  data.frame(xmin = boxes$x[i] - W/2 + 0.05 * k, xmax = boxes$x[i] + W/2 + 0.05 * k,
             ymin = boxes$y[i] - H/2 - 0.05 * k, ymax = boxes$y[i] + H/2 - 0.05 * k)))))
bell <- function(cx, cy, w = 0.42, h = 0.30) { t <- seq(-3, 3, length.out = 60)
  data.frame(x = cx + t / 6 * w, y = cy + h * exp(-t^2 / 2), cx = cx) }
G <- do.call(rbind, lapply(which(boxes$est), function(i) bell(boxes$x[i] + W/2 - 0.32, boxes$y[i] + H/2 + 0.08)))
seg <- function(x, y, xend, yend) data.frame(x, y, xend, yend)
A <- rbind(seg(6.0, 6.35, 8.3, 6.35), seg(8.3, 6.35, 8.3, 4.85),          # cell states -> mRNA
           seg(8.3, 3.62, 8.3, 1.75), seg(8.3, 1.75, 7.95, 1.75),          # mRNA -> protein (clear of the stacked outlines)
           seg(5.9, 1.75, 4.35, 1.75),                                      # protein -> binding
           seg(2.3, 1.75, 1.7, 1.75), seg(1.7, 1.75, 1.7, 3.58),           # binding -> communication
           seg(1.7, 4.85, 1.7, 6.35), seg(1.7, 6.35, 3.95, 6.35))          # communication -> cell states
head <- c(FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
loop <- data.frame(x = c(4.1, 3.85, 3.85, 4.45, 4.45), y = c(6.15, 6.15, 7.05, 7.05, 6.85))
g <- ggplot() + coord_fixed(xlim = c(0.4, 9.6), ylim = c(-0.35, 7.25), expand = FALSE) + theme_void() +
  geom_polygon(data = B, aes(x, y, group = g, fill = col, alpha = a), colour = NA) + scale_fill_identity() + scale_alpha_identity() +
  geom_path(data = loop, aes(x, y), colour = PUR, linewidth = 0.6, arrow = arrow(length = unit(0.09, "in"), type = "closed", ends = "first")) +
  geom_rect(data = stack, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", colour = RED, linewidth = 0.3) +
  geom_rect(data = boxes, aes(xmin = x - W/2, xmax = x + W/2, ymin = y - H/2, ymax = y + H/2), fill = "white", colour = PUR, linewidth = 0.5) +
  geom_text(data = boxes, aes(x, y, label = lab), family = FONT, size = pt(11), colour = PUR, lineheight = 0.9) +
  geom_segment(data = A[!head, ], aes(x, y, xend = xend, yend = yend), colour = PUR, linewidth = 0.6) +
  geom_segment(data = A[head, ], aes(x, y, xend = xend, yend = yend), colour = PUR, linewidth = 0.6,
               arrow = arrow(length = unit(0.09, "in"), type = "closed")) +
  geom_area(data = G, aes(x, y), fill = NA) +
  geom_polygon(data = G, aes(x, y, group = cx), fill = "#C9A8F0", colour = PUR, linewidth = 0.3) +
  geom_label(aes(x = 5.6, y = 4.55, label = "Spatialomics\ndata collection"), family = FONT, size = pt(11), colour = ORA, fill = "white",
             label.r = unit(0.15, "in"), label.padding = unit(0.08, "in"), lineheight = 0.9) +
  geom_label(aes(x = 5.0, y = 3.15, label = "Spatial modeling\n(posterior for every field)"), family = FONT, size = pt(10), colour = TEA, fill = "white",
             label.r = unit(0.15, "in"), label.padding = unit(0.08, "in"), lineheight = 0.9) +
  geom_segment(aes(x = 3.2, y = 1.05, xend = 4.95, yend = 0.25), colour = TEA, linewidth = 0.6, arrow = arrow(length = unit(0.09, "in"), type = "closed")) +
  geom_text(aes(x = 3.95, y = 0.45, label = "joint posterior"), family = FONT, size = pt(8), colour = TEA, hjust = 1) +
  geom_label(aes(x = 6.3, y = 0.2, label = "Credible differences:\nregions, transects, whole fields"), family = FONT, size = pt(10), colour = TEA, fill = "white",
             label.r = unit(0.05, "in"), label.padding = unit(0.08, "in"), lineheight = 0.9) +
  geom_text(aes(x = 0.5, y = -0.25, label = "Stacked outlines: high-dimensional (many genes or ligand-receptor pairs).  Bell glyph: posterior distribution of an estimated quantity."),
            family = FONT, size = pt(8), colour = "grey30", hjust = 0)
agg_png("fig1-schematic.png", width = 7.5, height = 7.5 * 7.6 / 9.2, units = "in", res = 600); print(g); invisible(dev.off())
cat("ok\n")
