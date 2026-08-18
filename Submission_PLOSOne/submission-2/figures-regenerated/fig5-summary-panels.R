## Rebuild only Figure 5's two summary panels, on a dark ground to match the maps.
##
## regenerate-fig5.R produced these with theme_classic() + DarkTheme(), which left them on a
## white ground: theme_classic sets panel.background after DarkTheme has run, so the dark
## elements are overridden. They are the only light panels in an otherwise dark figure.
## Split out here because they need only the small Sup6 objects, not the five large Seurat
## assays that the map panels require.
##
## Run: distrobox enter emacs-r -- Rscript fig5-summary-panels.R

suppressMessages({library(ggplot2); library(tidyr); library(ggrepel)})

repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
sup6 <- file.path(repo, "Supplements for Publication/Sup6 - ELSA Analysis")
pdir <- file.path(repo, "Submission_PLOSOne/submission-2/figures-regenerated/panels-fig5")

dark <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.background  = element_rect(fill = "black", colour = NA),
          panel.background = element_rect(fill = "black", colour = NA),
          panel.grid.major = element_line(colour = "grey22", linewidth = 0.25),
          panel.grid.minor = element_blank(),
          axis.text  = element_text(colour = "grey85"),
          axis.title = element_text(colour = "grey90"),
          axis.line  = element_line(colour = "grey70", linewidth = 0.3),
          legend.background = element_rect(fill = "black", colour = NA),
          legend.key = element_rect(fill = "black", colour = NA),
          legend.text  = element_text(colour = "grey85"),
          legend.title = element_text(colour = "grey90"))
}

e1 <- new.env(); load(file.path(sup6, "assay.mean.elsa.trends.2025-07-27.Robj"), envir = e1)
trends <- get(ls(e1)[1], envir = e1)
e2 <- new.env(); load(file.path(sup6, "scatter.dat.2025-08-07.Robj"), envir = e2)
scat <- as.data.frame(get(ls(e2)[1], envir = e2))

keep <- c("l.est", "r.est", "m.est.min", "m.est.prod", "m.est.gm", "l", "r")
td <- data.frame(trends); td$param_d <- seq(50, 1000, 50)
td <- pivot_longer(td, cols = -param_d, names_to = "assay", values_to = "elsa")
td <- td[td$assay %in% keep, ]

p_d <- ggplot(td, aes(param_d, elsa, colour = assay, group = assay)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.1) + ylim(0, NA) +
  labs(x = "ELSA neighbourhood radius (param_d)", y = "mean ELSA", colour = NULL) +
  scale_colour_brewer(palette = "Set2") + dark()
ggsave(file.path(pdir, "Fig5_summary_elsa-vs-paramD.png"), p_d,
       width = 6.5, height = 3.4, dpi = 300, bg = "black")

p_n <- ggplot(scat, aes(num.clust, elsa)) +
  geom_point(size = 2.4, colour = "grey70") +
  geom_point(data = scat[scat$assay == "m.est.min", ], size = 3.4, colour = "#E4572E") +
  geom_text_repel(data = scat[scat$assay %in% keep, ], aes(label = assay),
                  size = 2.5, colour = "grey85", max.overlaps = 20, segment.colour = "grey45") +
  ylim(0, NA) + labs(x = "number of clusters", y = "mean ELSA (d = 100)") + dark()
ggsave(file.path(pdir, "Fig5_summary_elsa-vs-numclust.png"), p_n,
       width = 6.5, height = 3.4, dpi = 300, bg = "black")

cat("rebuilt both Figure 5 summary panels on a dark ground\n")
