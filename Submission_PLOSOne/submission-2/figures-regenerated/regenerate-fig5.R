## Regenerate the Figure 5 panels (convolution operator comparison) at a correct 1:1
## aspect ratio.
##
## Why (R2.6, extended): Reviewer 2 flagged elongated pixels in Figs 3 and 4, caused by
## `asp = qp.res.y/qp.res.x` in base-graphics calls. Fig 5 has the same defect from a
## different source: its spatial panels come from ggplot in `Sup6 - ELSA Analysis.R`, which
## never calls coord_fixed(), so ggplot stretched each panel to fill whatever device was
## requested (7 x 5 in) rather than to the true 200:130 bin geometry. Adding coord_fixed()
## fixes it. The reviewer did not flag this one; the authors chose to correct it anyway.
##
## Nothing is re-clustered. Every panel is drawn from the stored cluster labels and stored
## fields of the 2025-07-25/26 run, so the published result is preserved exactly.
##
## Revised 2026-10-08: panels are rendered at the size they print in the transposed composite
## (compose-fig-5-transposed.py, 6.69 in wide), so type is 8 pt in print as PLOS requires
## (it was about 3 pt after shrinking). Titles that only repeated the composite's row and
## column labels are dropped. The per-operator copies of the input panels are no longer
## written; the composite shows each input once (fig5-input-panels.R).
##
## Run: distrobox enter emacs-r -- Rscript regenerate-fig5.R

suppressMessages({
  library(Seurat); library(ggplot2); library(data.table)
  library(elsa); library(raster); library(cowplot); library(viridisLite)
})

repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
hd   <- file.path(repo, "sam_sandbox/HD Embeddings")
sup6 <- file.path(repo, "Supplements for Publication/Sup6 - ELSA Analysis")
odir <- file.path(repo, "Submission_PLOSOne/submission-2/figures-regenerated")
pdir <- file.path(odir, "panels-fig5")
dir.create(pdir, showWarnings = FALSE, recursive = TRUE)

OPERATORS <- c(prod = "m.est.prod", gm = "m.est.gm", min = "m.est.min")
MECHS <- list(c("WNT4","FZD6"), c("RSPO1","LGR6"), c("DHH","PTCH1"), c("VEGFA","KDR"))

## Published palette, from Sup6 / Sup7.
cols.epi   <- c('firebrick','steelblue','springgreen','purple','salmon','skyblue','navyblue','violetred')
cols.endo  <- c('#40376E','violetred','grey20','tomato','sandybrown','saddlebrown')
cols.mes   <- c('royalblue','plum4','lightgoldenrod','lawngreen','forestgreen','dimgray','deeppink','red2','paleturquoise1','palevioletred')
cols.myelo <- c('orchid4','purple4','plum1','olivedrab2','slateblue','mediumvioletred','sienna','orange','seagreen','lightseagreen','mediumpurple4','#F8F4A6')
cols.lymph <- c('#E0CA3C','#989FCE','#261132','#F37748','#AFFC41','#347FC4','#0B5351','#00A9A5',
                '#4E8098','#D81159','#AB3428','#01BAEF','#F42C04','#610345','#044B7F')
cols.global <- unique(c(cols.epi, cols.endo, cols.mes, cols.lymph, cols.myelo,
                        'white','yellow','#A88FAC','#D2D68D','#CBD4C2','red','blue','purple','magenta'))

load_assay <- function(a) { e <- new.env(); load(file.path(hd, paste0(a, "assay.2025-07-25.Robj")), envir = e); e$temp }

## The fix. Every spatial panel goes through this.
TXT <- 8
spatial_theme <- function(p, title = NULL) {
  p <- p + coord_fixed(ratio = 1)                    # <- the R2.6 correction for Fig 5
  if (!is.null(title)) p <- p + ggtitle(title)
  p + Seurat::DarkTheme() + Seurat::NoAxes() +
    theme(text = element_text(family = "Liberation Sans"),
          plot.background = element_rect(fill = "black", colour = NA),
          legend.background = element_rect(fill = "black", colour = NA),
          legend.title = element_blank(),
          plot.title = element_text(size = TXT, hjust = 0.5),
          legend.text = element_text(size = TXT),
          legend.key.width = unit(0.06, "in"), legend.key.height = unit(0.16, "in"),
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.03, "in"),
          plot.margin = margin(1, 1, 1, 1),
          panel.grid = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank())   # DarkTheme draws one; drop it so every panel matches
}

## final printed sizes in the transposed composite: interaction 1.48 in wide, maps 1.98 in
save_panel <- function(p, name, w = 1.98, h = 1.36) {
  ggsave(file.path(pdir, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "black")
  invisible(name)
}

field_panel <- function(meta, values, title = NULL, thresh = 0.25, option = "A") {
  d <- cbind(meta[, c("x","y")], v = as.numeric(values))
  spatial_theme(
    ggplot(d, aes(x, y, colour = v)) + geom_point(shape = 15, size = 0.06) +
      scale_colour_viridis_c(option = option, n.breaks = 3,
                             limits = c(0, thresh * max(d$v, na.rm = TRUE)),
                             oob = scales::squish),
    title)
}

message("loading ligand and receptor assays ...")
lig <- load_assay("l.est"); rec <- load_assay("r.est")
lig_counts <- lig@assays[["l.est"]]$counts; rownames(lig_counts) <- rownames(lig)
rec_counts <- rec@assays[["r.est"]]$counts; rownames(rec_counts) <- rownames(rec)
meta_lr <- lig@meta.data

## ---- top block: 4 mechanisms x 3 operators, each ligand + receptor + interaction --------
for (op in names(OPERATORS)) {
  a <- OPERATORS[[op]]
  message("operator ", a, " ...")
  obj <- load_assay(a)
  cnt <- obj@assays[[a]]$counts; rownames(cnt) <- rownames(obj); colnames(cnt) <- colnames(obj)
  meta <- obj@meta.data

  for (m in MECHS) {
    L <- m[1]; R <- m[2]; mech <- paste(L, R, sep = "-")
    if (!(mech %in% rownames(cnt))) { message("  skip ", mech, " (absent)"); next }
    thr <- if (mech == "WNT4-FZD6") 0.2 else 0.5
    save_panel(field_panel(meta, cnt[mech, ], NULL, thr),
               sprintf("Fig5_%s_%s_interaction", op, mech), w = 1.48, h = 0.86)
  }

  ## ---- middle: the tissue domains this operator produces -------------------------------
  dom <- data.frame(x = meta$x, y = meta$y, cluster = factor(meta$seurat_clusters))
  save_panel(spatial_theme(
      ggplot(dom, aes(x, y, colour = cluster)) + geom_point(shape = 15, size = 0.08) +
        scale_colour_manual(values = rep(cols.global, length.out = nlevels(dom$cluster))) +
        guides(colour = "none"),
      sprintf("%d domains", nlevels(dom$cluster))),
    sprintf("Fig5_%s_domains", op))

  ## ---- below that: the corresponding entrogram -----------------------------------------
  r  <- rasterFromXYZ(data.frame(x = meta$x, y = meta$y, z = as.numeric(dom$cluster)))
  ev <- elsa(r, d = 50, categorical = TRUE)
  ed <- as.data.frame(ev, xy = TRUE); names(ed)[3] <- "ELSA"
  save_panel(spatial_theme(
      ggplot(ed, aes(x, y, fill = ELSA)) + geom_raster() +
        scale_fill_viridis_c(option = "A", limits = c(0, 1), breaks = c(0, 0.5, 1)) +
        geom_point(alpha = 0)),
    sprintf("Fig5_%s_entrogram", op))

  rm(obj, cnt); gc(verbose = FALSE)
}

## Summary panels: see fig5-summary-panels.R, which owns them.

message("wrote ", length(list.files(pdir, pattern = "[.]png$")), " panels to ", pdir)
