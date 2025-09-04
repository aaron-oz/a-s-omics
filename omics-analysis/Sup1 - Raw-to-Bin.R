# 2024-dec-6
# aoz
# 2025-04-24 msbr mod
# 2025-07-27 MSBR mod
options(future.globals.maxSize = 10.0 * 1e9)

# This script loads in the raw raw raw counts, subsets to tissue boundary, aggregates to 10, 25, and 50 unit squares, and saves the output for later modeling.

require(glue)
require(ggplot2)
require(reshape2)
require(Matrix)
require(Seurat)
require(dplyr)
require(data.table)
require(splancs)
require(fmesher)
library(BPCells)

setwd("/Users/msbr/Dropbox/a-s-omics/data-inputs/mouse-embryo-raw/")

run.entire.preprocess <- FALSE

##############################################################
#### Load data (raw counts, unbinned) -- VERY LARGE FILE ~
##############################################################
## this is a data frame with columns: gene, x, y, and value
#if(run.entire.preprocess){
  full.raw <- read.delim("./full-raw-t16-data.tsv", comment.char="#") |> as.data.table()
  setnames(full.raw, 'geneID', 'feat')
  
  
  ## # Number of spots?
  ## nrow(E16.5_E1S3_GEM_bin1) # 552290433 measurements
  ## # Range of X and Y?
  ## range(E16.5_E1S3_GEM_bin1$x) # 0 44099
  ## range(E16.5_E1S3_GEM_bin1$y) # 0 26459
  
  # Limit to just ligand and receptor genes # msbr removing to act on all genes
  #lr.names <- readRDS('ligand-receptor-names.RDs') # msbr removing to act on all genes
  setkey(full.raw, feat)
  #sub <- subset(full.raw, feat %in% lr.names) # msbr removing to act on all genes
  sub <- full.raw
  nrow(full.raw) # 28044528
#}

# define boundary (found visually, see code below)
boundary <- matrix(ncol = 2, data = c(20000, 9500,
                                      19000, 9500,
                                      17000, 9600,
                                      13250, 10000,
                                      10500, 11000,
                                      8200, 13000,
                                      6250, 16000,
                                      6000, 18500,
                                      6800, 20800,
                                      7300, 22300,
                                      8700, 23400,
                                      10500, 24400,
                                      16000, 25500,
                                      19120, 22800,
                                      18900, 19200,
                                      20400, 19100,
                                      22000, 22000,
                                      24200, 22200,
                                      27800, 21800,
                                      29300, 22800,
                                      30600, 22600,
                                      31700, 21500,
                                      32200, 20800,
                                      32800, 19300,
                                      32800, 17500,
                                      32300, 15300,
                                      31400, 13700,
                                      29600, 11700,
                                      28500, 10800,
                                      27000, 9800,
                                      26500, 9600,
                                      25500, 9500,
                                      20000, 9500), byrow = T)

setwd("/Users/msbr/Dropbox/a-s-omics/data-inputs/mouse-embryo-raw/msbr-2025-04-24")

#if(run.entire.preprocess){
  feat.in.d <- which(splancs::inout(splancs::as.points(sub[, .(x, y)]), boundary))
  sub.bounded <- sub[feat.in.d, ]
  dim(sub.bounded) # 23660265, 15% reduction
  fwrite(sub.bounded, file.path('partial-raw-t16-data_2025-04-24.csv'))
  
  total.cts <- full.raw[, list(total.count = sum(MIDCount)), by = .(x, y)]
  total.in.d <- which(splancs::inout(splancs::as.points(total.cts[, .(x, y)]), boundary))
  total.cts <- total.cts[total.in.d, ]
  fwrite(total.cts, file.path('total-raw-t16-data_2025-04-24.csv'))
#}

## binning the data
gc()
# for(agg.size in c(25, 50, 75, 100)){
#   gc()
agg.size <- 50 # tester
  ## load in clean data
  sub.bounded <- fread('./partial-raw-t16-data_2025-04-24.csv')
  total.cts <- fread(file.path('total-raw-t16-data_2025-04-24.csv'))
  
  ## define aggregation bins
  
  x.breaks = seq(6000, 32800, by = agg.size)
  x.break.val = data.table(xbin = 1:(length(x.breaks) - 1),
                           x = x.breaks[-length(x.breaks)] + agg.size / 2)
  y.breaks = seq(9500, 25500, by = agg.size)
  y.break.val = data.table(ybin = 1:(length(y.breaks) - 1),
                           y = y.breaks[-length(y.breaks)] + agg.size / 2)
  
  ## aggregate feature
  
  print(glue('aggregating feature to agg.size {agg.size}'))
  
  sub.bounded[, xbin := cut(x, breaks = x.breaks, labels = FALSE)]
  sub.bounded[, ybin := cut(y, breaks = y.breaks, labels = FALSE)]
  all.dat <- sub.bounded[, list(feat.count = sum(MIDCount)),
                         by = .(xbin, ybin, feat)]
  # merge on x and y vals that are the center of the bins
  all.dat <- merge(all.dat, x.break.val,
                   by = 'xbin', all.x = T, all.y = F)
  all.dat <- merge(all.dat, y.break.val,
                   by = 'ybin', all.x = T, all.y = F)
  all.dat <- na.omit(all.dat)
  # again, drop point in case centers of bins are now outside boundary
  feat.in.d <- which(splancs::inout(splancs::as.points(all.dat[, .(x, y)]), boundary))
  all.dat <- all.dat[feat.in.d, ]
  
  # Save aggregated feature data
  fwrite(all.dat, file.path('all.dat.bin.50_2025-04-27.csv'))
  
  # Turn into DGE
  all.dat <- fread('./msbr-2025-04-24/all.dat.bin.50_2025-04-27.csv')
  # long to wide [had to break into parts to prvent data exceed]
  all.dat$coord <- paste(all.dat$x,all.dat$y,sep = '.') # add xy coordinate as a barcode
  split.by.gene <- split(all.dat,by='feat')
  feature.names <- names(split.by.gene)
  num.features <- length(feature.names)
  temp.dge.list <- list()
  for(i in 1:29){
    print(i)
    if(i!=29){
  to.run <- feature.names[((1000*i)-999):(1000*i)]
    }else{
      to.run <- feature.names[((1000*i)-999):num.features]
    }
    #temp.dat <- all.dat[all.dat$feat %in% to.run,]
    temp.dat <- subset(all.dat,subset = feat %in% to.run)
  temp.dge.list[[i]] <- dcast(temp.dat, coord ~ feat, value.var = "feat.count") # cast long to wide, using feature counts as the values
  temp.dge.list[[i]] <- as.data.frame(temp.dge.list[[i]])
  rownames(temp.dge.list[[i]]) <- temp.dge.list[[i]][,1] # put coord barcodes as rownames
  temp.dge.list[[i]] <- temp.dge.list[[i]][,-1] # remove barcodes from data table
  # Transpose for input into Seurat & convert to data.table for rbindlist below
  temp.dge.list[[i]] <- data.table(t(temp.dge.list[[i]]))
  }
  
  # Assemble into single DGE and save
  bin.50.dge <- rbindlist(temp.dge.list,fill = TRUE)
  rownames(bin.50.dge) <- feature.names # risky but i think OK #### DOESN'T DO ANYTHING b/c data.table does not allow rownames
  fwrite(bin.50.dge,file.path('bin.50.dge.csv'))
  
  # 2025-07-27 MSBR
  bin.50.dge <- fread('./msbr-2025-04-24/bin.50.dge.csv',header = T)
  bin.50.dge[1:5,1:5] # note that this does not have rownames, b/c is data.table
  gc()
  
  # Limit to features expressed in at least 0.25% of cells and at least 1 cell with two counts
  ncol(bin.50.dge) #126320
  
  foi <- feature.names[rowSums(bin.50.dge>0,na.rm = T) > (ncol(bin.50.dge)*0.0025) &
                         rowSums(bin.50.dge>=2,na.rm = T) > 0] # takes a long time and eats RAM
  length(foi) # 17120 (down from 28502)
  
  bin.50.dge$feature <- feature.names # adding feature column so that subsetting to foi can be performed, since no rownames..
  
  setkey(bin.50.dge, feature)
  
  temp <- subset(bin.50.dge, feature %in% foi)
  temp2 <- temp[,-c('feature')] # remove the feature column, so it is just a counts matrix
  temp3 <- Matrix::Matrix(as.matrix(temp2))
  
  # turn into seurat [cannot run because data exceeds 2^31 values, see https://github.com/satijalab/seurat/issues/9798]
  # bin.50.seurat <- CreateAssay5Object(counts = temp3)
  # bin.50.seurat <- CreateSeuratObject(bin.50.seurat,assay = 'stereoSeq')
  # temp@meta.data$x <- data.list[[1]]$x
  # temp@meta.data$y <- data.list[[1]]$y
  # temp
  
  # test <- CreateSeuratObject(counts = temp,assay = 'bin50')
  # 
  # # Save
  # bin.50.seurat <- test
  # save(bin.50.seurat,file = 'bin.50.seurat_2025-04-27.Robj')
  # 
  # 
  
  # But we really only need a few simple visualizations...
  nUMI <- colSums(temp2,na.rm = T)
  x.y <- colnames(temp2)
  
  meta.compiled <- data.frame(nUMI = nUMI,
                              x.y = x.y)
  meta.compiled$x <- stringr::str_split_fixed(meta.compiled$x.y,pattern = '[.]',n=2)[,1]
  meta.compiled$y <- stringr::str_split_fixed(meta.compiled$x.y,pattern = '[.]',n=2)[,2]
  
  fantom5 <- NICHES::ncomms8866_mousx.yfantom5 <- NICHES::ncomms8866_mouse
  ligand.list <- fantom5$Ligand.ApprovedSymbol
  receptor.list <- fantom5$Receptor.ApprovedSymbol
  ligand.rows <- which(temp$feature %in% ligand.list) # 452 features
  receptor.rows <- which(temp$feature %in% receptor.list) # 483 features
  
  temp.lig <- subset(temp, feature %in% ligand.list)
  temp.rec <- subset(temp, feature %in% receptor.list)
  
  temp.lig2 <- temp.lig[,-c('feature')] # remove the feature column, so it is just a counts matrix
  temp.rec2 <- temp.rec[,-c('feature')] # remove the feature column, so it is just a counts matrix

  nUMI.lig <- colSums(temp.lig2,na.rm = T)
  nUMI.rec <- colSums(temp.rec2,na.rm = T)

  meta.compiled$nUMI.lig <- nUMI.lig
  meta.compiled$nUMI.rec <- nUMI.rec
  
  temp.lig3 <- temp.lig2>0
  temp.rec3 <- temp.rec2>0
  
  nFeature.lig <- colSums(temp.lig3,na.rm = T)
  nFeature.rec <- colSums(temp.rec3,na.rm = T)
  
  meta.compiled$nFeature.lig <- nFeature.lig
  meta.compiled$nFeature.rec <- nFeature.rec
  
  meta.compiled$x <- as.numeric(meta.compiled$x)
  meta.compiled$y <- as.numeric(meta.compiled$y)
  
  save(meta.compiled,file = 'meta.compiled.MSBR.2025-07-27.Robj')

# plotting
setwd("~/Dropbox/a-s-omics/sam_sandbox/plotting for figures")

x.range <- c(15025,24975) # from endpoint seurat objects :: range(temp$x)
y.range <- c(11525,17975) # from endpoint seurat objects :: range(temp$y)

pdf(file = 'nUMI whole embryo.pdf',width=7, height = 4)
par(mar = c(0, 0, 0, 0))
ggplot(meta.compiled,
           aes(x=x,
               y=y,
               color = nUMI))+
      geom_point(size = 0.5,shape = 15)+
  annotate(geom = "rect", 
           ymax = y.range[2], 
           ymin = y.range[1], 
           xmax = x.range[2], 
           xmin = x.range[1], 
           colour = "white",
           linetype = 'dashed',
           fill = NA,
           linewidth = 0.75)+
      ylab(NULL)+
      xlab(NULL)+
      scale_color_viridis_c(option = 'E',NULL)+
      ggtitle('Total Counts')+
      DarkTheme()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
      NoAxes()
dev.off()


pdf(file = 'nUMI inset.pdf',width=7, height = 4)
par(mar = c(0, 0, 0, 0))
ggplot(meta.compiled[meta.compiled$x >= x.range[1] & meta.compiled$x <= x.range[2] &
                       meta.compiled$y >= y.range[1] & meta.compiled$y <= y.range[2],],
       aes(x=x,
           y=y,
           color = nUMI))+
  geom_point(size = 0.5,shape = 15)+
  ylab(NULL)+
  xlab(NULL)+
  scale_color_viridis_c(option = 'E',NULL)+
  ggtitle('Total Counts (inset)')+
  DarkTheme()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  NoAxes()
dev.off()

pdf(file = 'ligand nUMI.pdf',width=7, height = 4)
par(mar = c(0, 0, 0, 0))
ggplot(meta.compiled[meta.compiled$x >= x.range[1] & meta.compiled$x <= x.range[2] &
                       meta.compiled$y >= y.range[1] & meta.compiled$y <= y.range[2],],
       aes(x=x,
           y=y,
                        color = nUMI.lig))+
  geom_point(size = 0.5,shape = 15)+
  ylab(NULL)+
  xlab(NULL)+
  scale_color_viridis_c(option = 'B',NULL)+
  ggtitle('Total Ligand Counts')+
  DarkTheme()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  NoAxes()
dev.off()

pdf(file = 'receptor nUMI.pdf',width=7, height = 4)
par(mar = c(0, 0, 0, 0))
ggplot(meta.compiled[meta.compiled$x >= x.range[1] & meta.compiled$x <= x.range[2] &
                       meta.compiled$y >= y.range[1] & meta.compiled$y <= y.range[2],],
       aes(x=x,
           y=y,
                               color = nUMI.rec))+
  geom_point(size = 0.5,shape = 15)+
  ylab(NULL)+
  xlab(NULL)+
  scale_color_viridis_c(option = 'B',NULL)+
  ggtitle('Total Receptor Counts')+
  DarkTheme()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  NoAxes()
dev.off()

pdf(file = 'ligand nFeature.pdf',width=7, height = 4)
par(mar = c(0, 0, 0, 0))
ggplot(meta.compiled[meta.compiled$x >= x.range[1] & meta.compiled$x <= x.range[2] &
                       meta.compiled$y >= y.range[1] & meta.compiled$y <= y.range[2],],
       aes(x=x,
           y=y,
                               color = nFeature.lig))+
  geom_point(size = 0.5,shape = 15)+
  ylab(NULL)+
  xlab(NULL)+
  scale_color_viridis_c(option = 'D',NULL)+
  ggtitle('Number of Ligand Features Expressed')+
  DarkTheme()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  NoAxes()
dev.off()

pdf(file = 'receptor nFeature.pdf',width=7, height = 4)
par(mar = c(0, 0, 0, 0))
ggplot(meta.compiled[meta.compiled$x >= x.range[1] & meta.compiled$x <= x.range[2] &
                       meta.compiled$y >= y.range[1] & meta.compiled$y <= y.range[2],],
       aes(x=x,
           y=y,
                                 color = nFeature.rec))+
  geom_point(size = 0.5,shape = 15)+
  ylab(NULL)+
  xlab(NULL)+
  scale_color_viridis_c(option = 'D',NULL)+
  ggtitle('Number of Receptor Features Expressed')+
  DarkTheme()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  NoAxes()
dev.off()

pdf(file = 'ligand nUMI fraction.pdf',width=7, height = 4)
par(mar = c(0, 0, 0, 0))
ggplot(meta.compiled[meta.compiled$x >= x.range[1] & meta.compiled$x <= x.range[2] &
                       meta.compiled$y >= y.range[1] & meta.compiled$y <= y.range[2],],
       aes(x=x,
           y=y,
                                   color = nUMI.lig/nUMI))+
  geom_point(size = 0.5,shape = 15)+
  ylab(NULL)+
  xlab(NULL)+
  scale_color_viridis_c(option = 'G',NULL)+
  ggtitle('Fraction of counts that are ligands')+
  DarkTheme()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  NoAxes()
dev.off()

pdf(file = 'receptor nUMI fraction.pdf',width=7, height = 4)
par(mar = c(0, 0, 0, 0))
ggplot(meta.compiled[meta.compiled$x >= x.range[1] & meta.compiled$x <= x.range[2] &
                       meta.compiled$y >= y.range[1] & meta.compiled$y <= y.range[2],],
       aes(x=x,
           y=y,
           color = nUMI.rec/nUMI))+
  geom_point(size = 0.5,shape = 15)+
  ylab(NULL)+
  xlab(NULL)+
  scale_color_viridis_c(option = 'G',NULL)+
  ggtitle('Fraction of counts that are receptors')+
  DarkTheme()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  NoAxes()
dev.off()

pdf(file = 'receptor nUMI fraction.pdf',width=7, height = 4)
par(mar = c(0, 0, 0, 0))
ggplot(meta.compiled[meta.compiled$x >= x.range[1] & meta.compiled$x <= x.range[2] &
                       meta.compiled$y >= y.range[1] & meta.compiled$y <= y.range[2],],
       aes(x=x,
           y=y,
           color = nUMI.rec/nUMI))+
  geom_point(size = 0.5,shape = 15)+
  ylab(NULL)+
  xlab(NULL)+
  scale_color_viridis_c(option = 'G',NULL)+
  ggtitle('Fraction of counts that are receptors')+
  DarkTheme()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  NoAxes()
dev.off()

                        