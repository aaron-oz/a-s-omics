## Build the controls figure (S10 Fig) for reviewer comment R2.8.
##
## Panels
##   A  spatial permutation test, domain coherence: null distribution of mean ELSA(d50)
##      across 49 coordinate permutations, each a full refit of all 912 features
##   B  the same permutation test, number of domains recovered. "Domains" is the count of
##      distinct Louvain community labels returned for that permutation
##   C  the emergent signalling domains
##   D  the anatomical annotation published for the same section
##   E  per-domain agreement: purity of each emergent domain against its dominant label
##
## Colour: each anatomical label gets one hue. Every emergent domain in C is drawn in the
## hue of the label it best matches, with luminance varied when several domains map to the
## same label. Corresponding regions therefore share a hue between C and D, while domains
## that subdivide one organ stay distinguishable.
##
## Run: distrobox enter emacs-r -- Rscript make-controls-figure.R
suppressMessages(library(data.table))
repo <- "/var/home/aoz/Dropbox/genetics/a-s-omics"
cd   <- file.path(repo, "data-outputs/mouse-embryo-raw/controls")
odir <- file.path(repo, "Submission_PLOSOne/submission-2/figures-regenerated")

B   <- rbindlist(lapply(Sys.glob(file.path(cd,"control-B-klone","control-B-perm*.csv")), fread))
obs <- B[perm_id==0]; nul <- B[perm_id!=0]
C   <- readRDS(file.path(cd,"control-C-result.rds"))
lab <- readRDS(file.path(cd,"observed-domain-labels.rds")); setnames(lab,"domain","dom")
ann <- fread(file.path(cd,"e16-fov-annotation.csv"))
D   <- merge(ann, lab, by=c("x","y"))[!is.na(annotation)]

## ---- shared colour scheme ----
pur  <- as.data.table(C$purity); setorder(pur, -purity)
labs <- names(sort(table(D$annotation), decreasing=TRUE))
## golden-angle hue assignment: labels are ordered by frequency, so stepping 137.5 degrees
## puts the LARGEST regions on maximally separated hues instead of adjacent ones
hue  <- setNames((15 + 137.508*(seq_along(labs)-1)) %% 360, labs)
ann.hex <- setNames(hcl(h=hue[labs], c=72, l=62), labs)
dom.base <- setNames(as.character(pur$dominant), as.character(pur$domain))
dom.hex  <- character(0)
for (lb in unique(dom.base)) {
  ds <- names(dom.base)[dom.base==lb]
  ll <- if (length(ds)==1) 62 else seq(42, 80, length.out=length(ds))
  dom.hex[ds] <- hcl(h=hue[[lb]], c=72, l=ll)
}

## exact categorical rendering; quilt.plot bins z onto a ramp and mis-assigns codes
ux <- sort(unique(D$x)); uy <- sort(unique(D$y))
cellmap <- function(code){ M <- matrix(NA_real_, length(ux), length(uy))
  M[cbind(match(D$x,ux), match(D$y,uy))] <- code; M }
catplot <- function(M, cols, ttl){
  image(ux, uy, M, col=cols, breaks=seq(0.5,length(cols)+0.5,1),
        main=ttl, xlab="", ylab="", useRaster=TRUE, asp=1); box() }

png(file.path(odir,"S10_Fig_controls.png"), width=13, height=11.5, units="in", res=300)
layout(matrix(c(1,2, 3,4, 5,5), nrow=3, byrow=TRUE), heights=c(1, 1.45, 1))

## A ----
par(mar=c(4.2,4.4,3.4,1.4), cex.lab=1.05)
h <- hist(nul$elsa_d50, breaks=12, plot=FALSE)
plot(h, col="grey80", border="white", xlim=c(obs$elsa_d50*0.5, max(nul$elsa_d50)*1.04),
     main="A   Spatial permutation test: domain coherence",
     xlab="mean ELSA at radius 50", ylab="permutations")
abline(v=obs$elsa_d50, col="firebrick", lwd=3)
text(obs$elsa_d50, max(h$counts)*0.92, sprintf("  observed %.3f", obs$elsa_d50), col="firebrick", adj=0, font=2)
text(mean(nul$elsa_d50), max(h$counts)*0.55, sprintf("49 permutations\nmean %.3f", mean(nul$elsa_d50)), adj=0.5, col="grey25")
mtext("lower = more spatially coherent", side=3, line=0.1, cex=0.72, col="grey35")

## B ----
h2 <- hist(nul$n_clusters, breaks=seq(min(nul$n_clusters)-0.5, max(nul$n_clusters)+0.5,1), plot=FALSE)
plot(h2, col="grey80", border="white", xlim=c(min(nul$n_clusters)-1, obs$n_clusters+2),
     main="B   Spatial permutation test: number of domains",
     xlab="domains recovered (distinct cluster labels)", ylab="permutations")
abline(v=obs$n_clusters, col="firebrick", lwd=3)
text(obs$n_clusters, max(h2$counts)*0.92, sprintf("  observed %d", obs$n_clusters), col="firebrick", adj=0, font=2)

## C and D ----
par(mar=c(3.2,3.2,3.4,0.8))
dlev <- as.character(pur$domain)
catplot(cellmap(match(as.character(D$dom), dlev)), unname(dom.hex[dlev]), "C   Emergent signalling domains")
catplot(cellmap(match(as.character(D$annotation), labs)), unname(ann.hex[labs]),
        "D   Published anatomical annotation, same section")

## E ----
par(mar=c(9.8,4.6,3.4,1.4))
bp <- barplot(100*pur$purity, col=unname(dom.hex[as.character(pur$domain)]), border="grey30",
              ylim=c(0,105), ylab="purity against dominant label (%)",
              main="E   Per-domain agreement with the anatomical annotation",
              names.arg=rep("",nrow(pur)))
text(bp, -4, srt=45, adj=1, xpd=NA, cex=0.78, labels=sprintf("%s (n=%d)", pur$dominant, pur$n))
wm <- 100*sum(pur$purity*pur$n)/sum(pur$n)
abline(h=wm, lty=2, col="firebrick", lwd=2)
text(bp[length(bp)], wm+5, sprintf("size-weighted mean %.1f%%", wm), col="firebrick", adj=1, font=2, cex=0.9)
mtext(sprintf("overall ARI %.3f, NMI %.3f, permutation p = %.3f (499 permutations, all ARI < %.4f)",
              C$ari, C$nmi, C$p, max(C$null)), side=3, line=0.1, cex=0.74, col="grey35")
dev.off()
cat(sprintf("wrote %s\n", file.path(odir,"S10_Fig_controls.png")))
