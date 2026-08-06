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
spatial_theme <- function(p, title) {
  p + coord_fixed(ratio = 1) +                       # <- the R2.6 correction for Fig 5
    ggtitle(title) + Seurat::DarkTheme() + Seurat::NoAxes() +
    theme(legend.title = element_blank(),
          plot.title = element_text(size = 10),
          panel.grid = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank())   # DarkTheme draws one; drop it so every panel matches
}

save_panel <- function(p, name, w = 5.2, h = 3.6) {
  ggsave(file.path(pdir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "black")
  invisible(name)
}

field_panel <- function(meta, values, title, thresh = 0.25, option = "A") {
  d <- cbind(meta[, c("x","y")], v = as.numeric(values))
  spatial_theme(
    ggplot(d, aes(x, y, colour = v)) + geom_point(shape = 15, size = 0.45) +
      scale_colour_viridis_c(option = option,
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
    save_panel(field_panel(meta_lr, lig_counts[L, ], paste(L, "| l.est"), thr),
               sprintf("Fig5_%s_%s_ligand", op, mech), w = 3.0, h = 2.1)
    save_panel(field_panel(meta_lr, rec_counts[R, ], paste(R, "| r.est"), thr),
               sprintf("Fig5_%s_%s_receptor", op, mech), w = 3.0, h = 2.1)
    save_panel(field_panel(meta, cnt[mech, ], paste(mech, "|", a), thr),
               sprintf("Fig5_%s_%s_interaction", op, mech), w = 5.2, h = 3.6)
  }

  ## ---- middle: the tissue domains this operator produces -------------------------------
  dom <- data.frame(x = meta$x, y = meta$y, cluster = factor(meta$seurat_clusters))
  save_panel(spatial_theme(
      ggplot(dom, aes(x, y, colour = cluster)) + geom_point(shape = 15, size = 0.45) +
        scale_colour_manual(values = rep(cols.global, length.out = nlevels(dom$cluster))) +
        guides(colour = "none"),
      sprintf("%s | %d domains", a, nlevels(dom$cluster))),
    sprintf("Fig5_%s_domains", op))

  ## ---- below that: the corresponding entrogram -----------------------------------------
  r  <- rasterFromXYZ(data.frame(x = meta$x, y = meta$y, z = as.numeric(dom$cluster)))
  ev <- elsa(r, d = 50, categorical = TRUE)
  ed <- as.data.frame(ev, xy = TRUE); names(ed)[3] <- "ELSA"
  save_panel(spatial_theme(
      ggplot(ed, aes(x, y, fill = ELSA)) + geom_raster() +
        scale_fill_viridis_c(option = "A", limits = c(0, 1)) +
        geom_point(alpha = 0),
      sprintf("%s | entrogram, ELSA d=50", a)),
    sprintf("Fig5_%s_entrogram", op))

  rm(obj, cnt); gc(verbose = FALSE)
}

## ---- bottom row: the two summary panels, no spatial aspect to correct -------------------
message("summary panels ...")
e1 <- new.env(); load(file.path(sup6, "assay.mean.elsa.trends.2025-07-27.Robj"), envir = e1)
trends <- get(ls(e1)[1], envir = e1)
e2 <- new.env(); load(file.path(sup6, "scatter.dat.2025-08-07.Robj"), envir = e2)
scat <- get(ls(e2)[1], envir = e2)

keep <- c("l.est","r.est","m.est.min","m.est.prod","m.est.gm","l","r")
td <- data.frame(trends); td$param_d <- seq(50, 1000, 50)
td <- tidyr::pivot_longer(td, cols = -param_d, names_to = "assay", values_to = "elsa")
td <- td[td$assay %in% keep, ]

p_d <- ggplot(td, aes(param_d, elsa, colour = assay, group = assay)) +
  geom_point(size = 1.2) + geom_line() + ylim(0, NA) +
  labs(x = "ELSA neighbourhood radius (param_d)", y = "mean ELSA") +
  theme_classic() + Seurat::DarkTheme()
ggsave(file.path(pdir, "Fig5_summary_elsa-vs-paramD.png"), p_d,
       width = 6.5, height = 3.4, dpi = 300, bg = "black")

sc <- as.data.frame(scat)
p_n <- ggplot(sc, aes(num.clust, elsa)) +
  geom_point(size = 2.6, colour = "grey70") +
  geom_point(data = sc[sc$assay == "m.est.min", ], size = 3.4, colour = "red") +
  ggrepel::geom_text_repel(data = sc[sc$assay %in% keep, ], aes(label = assay),
                           size = 2.6, colour = "white", max.overlaps = 20) +
  ylim(0, NA) + labs(x = "number of clusters", y = "mean ELSA (d = 100)") +
  theme_classic() + Seurat::DarkTheme()
ggsave(file.path(pdir, "Fig5_summary_elsa-vs-numclust.png"), p_n,
       width = 6.5, height = 3.4, dpi = 300, bg = "black")

message("wrote ", length(list.files(pdir, pattern = "[.]png$")), " panels to ", pdir)
