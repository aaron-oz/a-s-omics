## Rebuild Figure 5's eight shared input panels (ligand and receptor fields).
##
## These sit in a narrow left gutter now that they are shown once rather than repeated in each
## operator column, and at that size an individual colour bar consumes roughly 40% of the panel
## while being unreadable. They are contextual, showing where each partner is expressed, not
## quantitative reads, so the bar is dropped and the map fills the gutter. The interaction
## fields beside them keep their bars.
##
## Loads only the l.est and r.est assays, not the three operator assays that the interaction
## panels need.
##
## Run: distrobox enter emacs-r -- Rscript fig5-input-panels.R

suppressMessages({library(Seurat); library(ggplot2)})

repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
hd   <- file.path(repo, "sam_sandbox/HD Embeddings")
pdir <- file.path(repo, "Submission_PLOSOne/submission-2/figures-regenerated/panels-fig5")

MECHS <- list(c("WNT4","FZD6"), c("RSPO1","LGR6"), c("DHH","PTCH1"), c("VEGFA","KDR"))

load_assay <- function(a) { e <- new.env(); load(file.path(hd, paste0(a, "assay.2025-07-25.Robj")), envir = e); e$temp }
lig <- load_assay("l.est"); rec <- load_assay("r.est")
lc <- lig@assays[["l.est"]]$counts; rownames(lc) <- rownames(lig)
rc <- rec@assays[["r.est"]]$counts; rownames(rc) <- rownames(rec)
meta <- lig@meta.data

panel <- function(values, title, thresh) {
  d <- cbind(meta[, c("x","y")], v = as.numeric(values))
  ggplot(d, aes(x, y, colour = v)) +
    geom_point(shape = 15, size = 0.45) +
    scale_colour_viridis_c(option = "A", guide = "none",
                           limits = c(0, thresh * max(d$v, na.rm = TRUE)),
                           oob = scales::squish) +
    coord_fixed(ratio = 1) + ggtitle(title) +
    Seurat::DarkTheme() + Seurat::NoAxes() +
    theme(plot.title = element_text(size = 9, colour = "white"),
          plot.margin = margin(1, 1, 1, 1),
          panel.grid = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.line = element_blank())
}

n <- 0
for (m in MECHS) {
  L <- m[1]; R <- m[2]; mech <- paste(L, R, sep = "-")
  thr <- if (mech == "WNT4-FZD6") 0.2 else 0.5
  ## one file per mechanism, not one per operator column: these panels are the reason the
  ## duplication was removed from the figure, so duplicating them on disk would be perverse
  ggsave(file.path(pdir, sprintf("Fig5_input_%s_ligand.png", mech)),
         panel(lc[L, ], L, thr), width = 2.4, height = 1.68, dpi = 300, bg = "black")
  ggsave(file.path(pdir, sprintf("Fig5_input_%s_receptor.png", mech)),
         panel(rc[R, ], R, thr), width = 2.4, height = 1.68, dpi = 300, bg = "black")
  n <- n + 2
}
cat("rewrote", n, "input panels without colour bars\n")
