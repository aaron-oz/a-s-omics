# Gene list for AOZ 2025-04-04

# pack
require(data.table)
require(dplyr)
require(fmesher)
require(reshape2)
require(ggplot2)
require(ggthemes)

# Set WD
setwd("/Users/msbr/Dropbox/a-s-omics/sam_sandbox")

# Load data (binned 50 E16.5, same as AOZ is using)
load("~/Dropbox/a-s-omics/data-inputs/mouse-embryo-raw/pre-fit-obj-binned-50.Rdata")

# Inspect data
class(all.dat)
plot(mesh)
plot(boundary.coords)
plot(coords,pch = '.')
plot(domain)
plot(xy.obs,pch = '.')

# load LR reference AOZ used - see AOZ script '0000-process-raw-data.R'
readRDS("~/Dropbox/a-s-omics/data-inputs/mouse-embryo-raw/ligand-receptor-names.Rds")

# check that it is sam's list [it is! only 974 genes - OK]
load("~/Large Files/chen_stereo_mouse_embryo/ligands.only.for.aoz.2024-10-27.Robj")
load("~/Large Files/chen_stereo_mouse_embryo/receptors.only.for.aoz.2024-10-27.Robj")
lr.list <- c(rownames(ligands.only),rownames(receptors.only)) # 974 genes total
unique(all.dat$feat) # 974 genes - OK

# bounding box edges from AOZ 2025-04-04
sub.bound <- splancs::as.points(matrix(c(15000, 25000, 25000, 15000, 15000,
                                         11500, 11500, 18000, 18000, 11500), ncol = 2))

in.sb <- which(splancs::inout(all.dat[, .(x, y)], sub.bound))
sub.dat <- all.dat[in.sb, ]
boundary.coords <- sub.bound
length(unique(sub.dat$feat)) # 973 genes OK

# long to wide
sub.dat$coord <- paste(sub.dat$x,sub.dat$y,sep = '.') # add xy coordinate as a barcode
dge <- dcast(sub.dat, coord ~ feat, value.var = "feat.count") # cast long to wide, using feature counts as the values
dim(dge)
rownames(dge) <- dge[,1] # but coord barcodes as rownames
dge <- dge[,-1] # remove barcodes from data table
dim(dge)

# Transpose for input into Seurat
counts <- t(dge)
dim(counts) # 973 26000 OK

# Compute percent expressivity for all genes
row.sums <- rowSums(counts>0)
head(row.sums)
total.spots <- ncol(counts)
frac.pos <- (row.sums/total.spots)*100
frac.pos

# explore distributions
to.plot <- data.frame(frac.pos)
ggplot(to.plot,
       aes(x=frac.pos))+geom_density()

# maximums?
row.max <- apply(counts, 1, max, na.rm=TRUE)
to.plot.2 <- data.frame(row.max)
ggplot(to.plot.2,
       aes(x=row.max))+geom_density()

# combine and scatter?
to.plot.3 <- data.frame(frac.pos = to.plot$frac.pos,
                        row.max = to.plot.2$row.max)
to.plot.3$gene <- rownames(to.plot)

ggplot(to.plot.3,
       aes(row.max,frac.pos)) +
  geom_point()+
  scale_y_log10()+
  scale_x_log10()

# plot all genes, with values on plot
for(i in 1:length(to.plot.3$gene)){
  print(i)
  foi <- to.plot.3$gene[i]
  data.to.plot <- data.frame(sub.dat)
  data.to.plot <- data.to.plot[data.to.plot$feat == foi,]
  thresh.factor <- 0.75
  png(filename = paste(signif(to.plot.3[to.plot.3$gene == foi,]$row.max,2),'_',foi,'.png',sep=''),width = 10,height = 7,units = 'in',res = 300)
  print(ggplot(data.to.plot,
         aes(x = x,
             y = y,
             fill = feat.count,
             stroke = NA,
             alpha = feat.count))+geom_point(size=1.5,shape=22)+
    theme_classic()+
    Seurat::NoAxes()+
    theme(panel.background = element_rect(fill='black'),plot.background = element_rect(fill='black'))+
    ggtitle(foi)+theme(plot.title = element_text(color = "red"))+
    scale_fill_gradient2(low="lightblue", high="red",
                          limits = c(0, max(data.to.plot$feat.count)*thresh.factor),
                          midpoint = max(data.to.plot$feat.count)*thresh.factor*0.5, 
                          oob = scales::squish)+
    scale_alpha(range = c(0,1))+
    annotate("text", x=15500, y=18100, 
             label= paste('Percent.Pos = ',signif(to.plot.3[to.plot.3$gene == foi,]$frac.pos,2)),color = 'white',hjust = 0)+
    annotate("text", x=15500, y=18300, 
             label= paste('Max.Count = ',signif(to.plot.3[to.plot.3$gene == foi,]$row.max,2)),color = 'white',hjust = 0)
  )
  dev.off()
}

# Define thresholds for genes to consider
row.max.thresh <- 2
frac.pos.thresh <- 0.25
to.use <- to.plot.3[to.plot.3$frac.pos>=frac.pos.thresh &
                      to.plot.3$row.max>=row.max.thresh,]
save(to.use,file = 'features.to.use.2025-04-04.Robj') # 918 features

# how many mechanisms?
mechs <- NICHES::ncomms8866_mouse
mechs.to.use <- mechs[mechs$Ligand.ApprovedSymbol %in% to.use$gene &
        mechs$Receptor.ApprovedSymbol %in% to.use$gene,]
save(mechs.to.use,file = 'mechs.to.use.2025-04-04.Robj') # 1511 mechanisms

# Make into Seurat object
counts.to.use <- counts[to.use$gene,]
e16 <- Seurat::CreateSeuratObject(counts = counts.to.use,assay = 'RNA')
e16





