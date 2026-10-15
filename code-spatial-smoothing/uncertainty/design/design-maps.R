## Design-phase maps (Phase 2 of PIVOT-PLAN-2026-10-14.org): anatomy and total counts only.
## Deliberately shows NO gene or signaling field and NO posterior quantity, so the regions and
## slices are chosen from anatomy and tissue coverage alone.
## Bin index: i = (x - 15025)/50 + 1 in 1..200, j = (y - 11525)/50 + 1 in 1..130.
## Run: OMP_NUM_THREADS=1 distrobox enter emacs-r -- Rscript design-maps.R
suppressMessages({library(data.table); library(ggplot2)})
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
a <- fread(file.path(repo, "data-outputs/mouse-embryo-raw/controls/e16-fov-annotation.csv"))
a[, `:=`(i = (x - 15025)/50 + 1, j = (y - 11525)/50 + 1)]
p <- readRDS(file.path(repo, "data-outputs/mouse-embryo-raw/2025-04-10/prediction-objects/Cd44-pred-poi.rds"))$lambda
tc <- p[, .(x, y, total.count)]
a <- merge(a, tc, by = c("x", "y"))
big <- a[, .N, by = annotation][N >= 500, annotation]
a[, ann := ifelse(annotation %in% big, annotation, "other")]
g1 <- ggplot(a, aes(i, j, fill = ann)) + geom_raster() + coord_fixed() +
  scale_x_continuous(breaks = seq(0, 200, 10), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 130, 10), expand = c(0, 0)) +
  theme_bw(9) + theme(panel.grid.major = element_line(color = "grey30", linewidth = 0.2),
                      panel.ontop = TRUE, panel.background = element_blank()) +
  labs(title = "MOSTA annotation, E16.5_E1S3 field of view (bin index)", fill = NULL)
ggsave("design-annotation.png", g1, width = 11, height = 7, dpi = 130)
g2 <- ggplot(a, aes(i, j, fill = log10(total.count))) + geom_raster() + coord_fixed() +
  scale_fill_viridis_c() +
  scale_x_continuous(breaks = seq(0, 200, 10), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 130, 10), expand = c(0, 0)) +
  theme_bw(9) + theme(panel.grid.major = element_line(color = "grey30", linewidth = 0.2),
                      panel.ontop = TRUE, panel.background = element_blank()) +
  labs(title = "Total UMI count per bin (log10)")
ggsave("design-totalcount.png", g2, width = 11, height = 7, dpi = 130)
print(a[, .(N = .N, med_total = as.numeric(median(total.count)), q10_total = quantile(total.count, 0.1),
            imin = min(i), imax = max(i), jmin = min(j), jmax = max(j)), by = annotation][order(-N)])
cat("bins with total.count < 100:", a[total.count < 100, .N], "\n")
