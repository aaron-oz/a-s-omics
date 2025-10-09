## 2025-07-12
## MSBR

## This script loads ASOMICS data and performs an experiment in which the feature number is slowly increased and the object re-clustered so that we can measure the following outputs:
## Cluster number (resolution fixed)
## Mean spatial entropy of the clusters (ELSA, categorical, dist parameter @ 50)
## UMAP Embedding, Spatial Embedding, ELSA Spatial

# Takes about a day to run, for loop is x30 at 1-40+ minutes per loop, depending on feature depth

# Pack
require(Seurat)
require(ggplot2)
require(elsa)
require(tmap)
require(tidyr)
# Load data
# load("~/Dropbox/a-s-omics/sam_sandbox/assay.list.2025-05-16.Robj")
load("~/Dropbox/a-s-omics/sam_sandbox/HD Embeddings/assay.list.2025-07-26.Robj")
# Inspect
names(assay.list)

# colors
gc()
cols.epi <- c('firebrick','steelblue','springgreen','purple','salmon','skyblue','navyblue','violetred')

cols.endo <- c('#40376E','violetred','grey20','tomato','sandybrown','saddlebrown')

cols.mes <- c('royalblue','plum4','lightgoldenrod','lawngreen','forestgreen','dimgray','deeppink','red2','paleturquoise1','palevioletred')

cols.myeloid <- c('orchid4','purple4','plum1','olivedrab2','slateblue','mediumvioletred','sienna','orange','seagreen','lightseagreen','mediumpurple4','#F8F4A6')

cols.lymphoid <- c('#E0CA3C','#989FCE','#261132',
                   '#F37748','#AFFC41','#347FC4',
                   '#0B5351','#00A9A5','#4E8098',
                   '#D81159','#AB3428','#01BAEF',
                   '#F42C04','#610345','#044B7F')

cols.global <- c(cols.epi,cols.endo,cols.mes,cols.lymphoid,cols.myeloid,'white','yellow','#A88FAC','#D2D68D','#CBD4C2','red','blue','purple','magenta')

# grab assay of interest & clear space
obj <- assay.list[['m.est.min']]
rm(assay.list)
gc()

# organize for downstream (Seurat v5 bug fix)
assay_v3 <- CreateAssayObject(counts = obj[["m.est.min"]]$counts)
obj[["RNA"]] <- assay_v3
DefaultAssay(obj) <- 'RNA'
obj <- ScaleData(obj,features = rownames(obj))
rm(assay_v3)

# rank features based on global variance
feature.ranks <- data.frame(rank = obj@assays$m.est.min@meta.data$var.features.rank,
                            feature = rownames(obj@assays$m.est.min))
feature.ranks <- feature.ranks[order(feature.ranks$rank),]

# stash for safety
obj.stash <- obj

# for loop to iteratively increase feature number

feat.num <- c(10,
              seq(from = 50,
                  to = nrow(feature.ranks),
                  by = 50),
              nrow(feature.ranks))

# set up parallelization
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)


for(i in 1:length(feat.num)){
#foreach(i=18:length(feat.num)) %dopar% {

  # print i
  print(i)
  
  # set features to use
  foi <- feature.ranks[1:feat.num[i],]$feature
  
  # create UMAP embedding in this feature space, without any PCA
  obj <- RunUMAP(obj,
                  features = foi,
                  dims = NULL,
                  reduction.name = paste('UMAP',length(foi),'feature',sep=''))
  
  # create neighborhood graph in full feature space, again without any PCA
  obj <- FindNeighbors(obj,
                        dims = NULL,
                        features = foi,
                        graph.name = c(paste(length(foi),'featureKNN',sep=''),
                                       paste(length(foi),'featureSNN',sep='')))
  
  # find clusters based on the above feature-space neighborhood graph (SNN)
  obj <- FindClusters(obj,
                       resolution = 1, # fixed, for this specific experiment
                       graph.name = paste('X',length(foi),'featureSNN',sep=''),
                      cluster.name = paste(length(foi),'featureClusters',sep=''))
  
  # plot the clusters in spatial coordinates
  to.plot <- data.frame(cluster = obj$seurat_clusters,
                        x = obj$x,
                        y = obj$y)
  p1 <- ggplot(to.plot,
               aes(x=x,
                   y=y,
                   colour = cluster))+
    geom_point(size = 0.5)+
    theme_classic()+
    DarkTheme()+
    NoLegend()+
    scale_color_manual(values = cols.global)+
    ggtitle(paste(length(foi),'Features | Spatial | Resolution = 1'))
  
  # plot the clusters in UMAP space
  p2 <- DimPlot(obj,
                label=T,
                label.color = 'white',
                cols = cols.global,
                reduction = paste('UMAP',length(foi),'feature',sep=''))+
    theme_classic()+
    DarkTheme()+
    ggtitle(paste(length(foi),'Features | UMAP | Resolution = 1'))
  
  # save plot as png
  setwd("/Users/msbr/Dropbox/a-s-omics/sam_sandbox")
  png(filename = paste(length(foi),'Feature Clustering.png'),width = 18,height = 6,units = 'in',res = 300)
  print(p1|p2)
  dev.off()

  # save metadata
  meta.data <- obj@meta.data
  save(meta.data,
       file = 'meta.data.feature.depth.experiment.2025-07-12.Robj')
  
  # save embedding metadata
  embedding.data <- obj@reductions
  save(embedding.data,
       file = 'embedding.data.feature.depth.experiment.2025-07-12.Robj')
  
  # save neighborhood graph data
  graph.data <- obj@graphs
  save(graph.data,
       file = 'graph.data.feature.depth.experiment.2025-07-12.Robj')
  
  if(i %in% c(1,5,10,15,20,25,length(feat.num))){
    message(paste('Saving ',i,'th ','object to disk...',sep=''))
    save(obj,
         file = 'object.feature.depth.experiment.2025-07-12.Robj')
  }
  
}

# Remake the plots but with colors fixed (can be done from saved objects above)

# make color palette unique
cols.global.2 <- unique(cols.global)

# plot cluster plots again, with new color palette, and without re-calculating the clusters or embeddings, etc.
for(i in 1:length(feat.num)){
  
  # print i
  print(i)
  
  # set features to use
  foi <- feature.ranks[1:feat.num[i],]$feature
  
  # plot the clusters in spatial coordinates
  to.plot <- data.frame(cluster = obj[[paste(length(foi),'featureClusters',sep='')]][,1],
                        x = obj$x,
                        y = obj$y)
  p1 <- ggplot(to.plot,
               aes(x=x,
                   y=y,
                   colour = cluster))+
    geom_point(size = 0.5)+
    theme_classic()+
    DarkTheme()+
    NoLegend()+
    scale_color_manual(values = cols.global.2)+
    ggtitle(paste(length(foi),'Features | Spatial | Resolution = 1'))
  
  # plot the clusters in UMAP space
  Idents(obj) <- paste(length(foi),'featureClusters',sep='')
  p2 <- DimPlot(obj,
                label=T,
                label.color = 'white',
                cols = cols.global.2,
                reduction = paste('UMAP',length(foi),'feature',sep=''))+
    theme_classic()+
    DarkTheme()+
    ggtitle(paste(length(foi),'Features | UMAP | Resolution = 1'))
  
  # save plot as png
  setwd("/Users/msbr/Dropbox/a-s-omics/sam_sandbox")
  png(filename = paste(length(foi),'Feature Clustering (Color Corrected).png'),width = 18,height = 6,units = 'in',res = 300)
  print(p1|p2)
  dev.off()
}

# ELSA Analysis 2025-07-13 MSBR

# elsa test script from https://cran.r-project.org/web/packages/elsa/vignettes/elsa.html
file <- system.file('external/lc_example.grd',package='elsa') 
lc <- raster(file)
plot(lc,main='Land cover: a categorical map')
elc <- elsa(lc, d=2000,categorical = T)
cl <- colorRampPalette(c('darkblue','yellow','red','black'))(100) # specifying a color scheme
plot(elc, col=cl, main='ELSA')

# make a raster object from each seurat dataset metadata...https://gis.stackexchange.com/questions/404443/create-raster-from-dataframe-containing-characters-as-values-in-r

assay.raw.elsa.values <- list()
assay.mean.elsa.trends <- list()

for(i in 1:length(feat.num)){
  
  print(i)
  
  # set features to use
  foi <- feature.ranks[1:feat.num[i],]$feature #important for naming
  
  # isolate clustering of interest
  temp <- obj@meta.data

  temp$seurat_clusters <- obj[[paste(length(foi),'featureClusters',sep='')]][,1] # this puts the clustering of interest in the 'seurat_clusters' column
  temp$seurat_clusters <- as.factor(temp$seurat_clusters )
  temp$seurat_clusters_numeric <- as.numeric(temp$seurat_clusters)
  
  # convert to raster
  r <- rasterFromXYZ(temp[,c("x","y","seurat_clusters_numeric")])
  r[] <- factor(levels(temp$seurat_clusters)[r[]])
  
  # run elsa for a variety of distance parameters
  dist.list <-seq(50,1000,50)
  elsa.by.dist <- list()
  for(j in 1:length(dist.list)){
    print(paste('i = ',i,'/',length(feat.num), ' j = ',j,'/',length(dist.list),sep=''))
    test <- elsa(r, d=dist.list[j],categorical = T)
    elsa.by.dist[[j]] <- test
    
    assay.of.interest <- paste(length(foi),'featureClusters',sep='')
    
    setwd("/Users/msbr/Dropbox/a-s-omics/sam_sandbox")
    
    png(filename = paste(assay.of.interest,' ELSA | d = ',dist.list[j],' | Enterogram', '.png',sep = ''),
        width = 7,height = 5,res = 300,units = 'in')
    print(plot(elsa.by.dist[[j]], col=cl, main='ELSA'))
    dev.off()

    png(filename = paste(assay.of.interest,' ELSA | d = ',dist.list[j],' | Distribution', '.png',sep = ''),
        width = 7,height = 5,res = 300,units = 'in')
    print(ggplot(data = data.frame(elsa = elsa.by.dist[[j]]@data@values),
                 aes(x = elsa))+
            geom_density()+
            theme_classic()+
            ggtitle(paste('ELSA Distribution | d =',dist.list[[j]])))
    dev.off()

  }
  
  # Stash raw elsa values
  assay.raw.elsa.values[[i]] <- elsa.by.dist
  
  # Look at mean trend with d
  mean.elsa <- list()
  for(k in 1:length(elsa.by.dist)){
    mean.elsa[[k]] <- mean(elsa.by.dist[[k]]@data@values)
  }
  
  # Stash mean elsa values (trend with increasing d)
  assay.mean.elsa.trends[[i]] <- unlist(mean.elsa)
  assay.mean.elsa.trends[[i]] <- data.frame(elsa = assay.mean.elsa.trends[[i]],
                                            d = dist.list)
  # list wise naming
  names(assay.raw.elsa.values)[i] <- paste(length(foi),'featureClusters',sep='')
  names(assay.raw.elsa.values[[i]]) <- paste('d = ',dist.list,sep='')
  names(assay.mean.elsa.trends)[i] <- paste(length(foi),'featureClusters',sep='')
}

# save for later
save(assay.raw.elsa.values,file = 'clustering.by.nfeatures.raw.elsa.values.Robj')
save(assay.mean.elsa.trends,file = 'clustering.by.nfeatures.mean.elsa.trends.Robj')


# summary plotting
elsa.d50 <- c()
for(i in 1:length(feat.num)){
  elsa.d50 <- c(elsa.d50,assay.mean.elsa.trends[[i]][1,]$elsa)
}
elsa.d50 <- data.frame(feat.num = feat.num,
                       elsa.value = elsa.d50,
                       elsa.dist = '50')
elsa.d100 <- c()
for(i in 1:length(feat.num)){
  elsa.d100 <- c(elsa.d100,assay.mean.elsa.trends[[i]][2,]$elsa)
}
elsa.d100 <- data.frame(feat.num = feat.num,
                       elsa.value = elsa.d100,
                       elsa.dist = '100')
elsa.d150 <- c()
for(i in 1:length(feat.num)){
  elsa.d150 <- c(elsa.d150,assay.mean.elsa.trends[[i]][3,]$elsa)
}
elsa.d150 <- data.frame(feat.num = feat.num,
                        elsa.value = elsa.d150,
                        elsa.dist = '150')
elsa.d200 <- c()
for(i in 1:length(feat.num)){
  elsa.d200 <- c(elsa.d200,assay.mean.elsa.trends[[i]][4,]$elsa)
}
elsa.d200 <- data.frame(feat.num = feat.num,
                        elsa.value = elsa.d200,
                        elsa.dist = '200')

to.plot.trends <- rbind(elsa.d50,
                 elsa.d100,
                 elsa.d150)
to.plot.trends$elsa.dist <- factor(to.plot.trends$elsa.dist,
                               levels = c('50','100','150'))
ggplot(to.plot.trends,
       aes(x = feat.num,
           y = elsa.value,
           color = elsa.dist))+
  ylab('Spatial Entropy (ELSA)')+
  geom_point()+geom_line()+theme_classic()

# go get the variance score for each feature
feat.var.scores <- obj@assays$m.est.min@meta.data
ggplot(feat.var.scores,
       aes(x = vf_vst_counts_rank,
           y = vf_vst_counts_variance.standardized))+
  geom_bar(stat='identity')+
  theme_classic()

# go get cluster number for each feature set
num.clust <- data.frame()
for(i in 1:length(feat.num)){
  assay.of.interest <- paste(feat.num[i],'featureClusters',sep='')
  temp <- data.frame(num.of.features = feat.num[i],
                    num.of.clusters = length(table(obj[[assay.of.interest]])))
  num.clust <- rbind(num.clust,temp)
}
ggplot(num.clust,
       aes(x=num.of.features,
           y=num.of.clusters))+
geom_point()+geom_line()+theme_classic()

# assemble
p1 <- ggplot(to.plot,
             aes(x = feat.num,
                 y = elsa.value,
                 color = elsa.dist))+
  #ylim(0,0.3)+
  ylab('Spatial Entropy')+
  xlab('Features')+
  geom_point()+geom_line()+theme_classic()
p2 <- ggplot(feat.var.scores,
             aes(x = vf_vst_counts_rank,
                 y = vf_vst_counts_variance.standardized))+
  geom_bar(stat='identity')+
  ylab('Variance')+
  xlab('Features')+
  theme_classic()
p3 <- ggplot(num.clust,
             aes(x=num.of.features,
                 y=num.of.clusters))+
  geom_point()+
  geom_line()+
  ylab('Number of Clusters')+
  xlab('Features')+
  theme_classic()

# Save summary plot
setwd("/Users/msbr/Dropbox/a-s-omics/sam_sandbox")
png(filename = 'm.est.min feature depth experiment summary trends.png',
    width = 5,height = 8,res = 300,units = 'in')
print(cowplot::plot_grid(p1,p2,p3,ncol = 1,align = 'hv'))
dev.off()

# Final for loop, plotting spatial embedding, UMAP embedding, ELSA with d = 100, and the three trends all as a whole

for(i in 1:length(feat.num)){
  
  # print i
  print(i)
  
  # set features to use
  foi <- feature.ranks[1:feat.num[i],]$feature
  
  # plot the clusters in spatial coordinates
  to.plot.clusters <- data.frame(cluster = obj[[paste(length(foi),'featureClusters',sep='')]][,1],
                        x = obj$x,
                        y = obj$y)
  p1 <- ggplot(to.plot.clusters,
               aes(x=x,
                   y=y,
                   colour = cluster))+
    geom_point(size = 0.5)+
    theme_classic()+
    DarkTheme()+
    NoLegend()+
    scale_color_manual(values = cols.global.2)+
    ggtitle(paste(length(foi),'Features | Spatial | Resolution = 1'))
  
  # plot the clusters in UMAP space
  Idents(obj) <- paste(length(foi),'featureClusters',sep='')
  p2 <- DimPlot(obj,
                label=T,
                label.color = 'white',
                cols = cols.global.2,
                reduction = paste('UMAP',length(foi),'feature',sep=''))+
    theme_classic()+
    DarkTheme()+
    ggtitle(paste(length(foi),'Features | UMAP | Resolution = 1'))
  
  # elsa plot d = 100 https://tmieno2.github.io/R-as-GIS-for-Economists/geom-raster.html
  # first convert raster to df
  elsa.df <- as.data.frame(assay.raw.elsa.values[[i]][[1]],xy=T)#first index (i) = feat.num, second index (2) means elsa.d = 50
  p3 <- ggplot(data = elsa.df) +
    geom_raster(aes(x = x, 
                    y = y, 
                    fill = ELSA)) +
    scale_fill_viridis_c(option = 'A') +
    theme_void() +
    theme(legend.position = "bottom")
  p3

  # summary plot, with progression
  
  # progression of elsa
  p.a <- ggplot(to.plot.trends,
               aes(x = feat.num,
                   y = elsa.value,
                   color = elsa.dist))+
    ylab('Spatial Entropy')+
    xlab('Features')+
    geom_point()+
    geom_line()+
    theme_classic()+
    scale_color_manual(values = c('lightgrey','lightgrey','lightgrey'))+
    ggtitle('Spatial Coherence vs. Feature Number | Assay = m.est.min')
    p.a.i <- p.a +
    geom_point(data = to.plot.trends[to.plot.trends$feat.num<=feat.num[i] & to.plot.trends$elsa.dist == '50',],
               color = 'black')+
    geom_line(data = to.plot.trends[to.plot.trends$feat.num<=feat.num[i] & to.plot.trends$elsa.dist == '50',],
               color = 'black')+
    geom_point(data = to.plot.trends[to.plot.trends$feat.num<=feat.num[i] & to.plot.trends$elsa.dist == '100',],
               color = 'red')+
    geom_line(data = to.plot.trends[to.plot.trends$feat.num<=feat.num[i] & to.plot.trends$elsa.dist == '100',],
              color = 'red')+
    geom_point(data = to.plot.trends[to.plot.trends$feat.num<=feat.num[i] & to.plot.trends$elsa.dist == '150',],
               color = 'purple')+
    geom_line(data = to.plot.trends[to.plot.trends$feat.num<=feat.num[i] & to.plot.trends$elsa.dist == '150',],
              color = 'purple')
  
  # progression of features
  p.b <- ggplot(feat.var.scores,
               aes(x = vf_vst_counts_rank,
                   y = vf_vst_counts_variance.standardized))+
    geom_bar(stat='identity',color = 'lightgrey')+
    ylab('Variance')+
    xlab('Features')+
    theme_classic()

  p.b.i <- p.b + 
    geom_bar(data = feat.var.scores[feat.var.scores$var.features.rank<=feat.num[i],],
             stat = 'identity',
                          color = 'red')

  # progression of clusters
  p.c <- ggplot(num.clust,
               aes(x=num.of.features,
                   y=num.of.clusters))+
    geom_point(color = 'lightgrey')+
    geom_line(color = 'lightgrey')+
    ylab('Number of Clusters')+
    xlab('Features')+
    theme_classic()
  
  p.c.i <- p.c + 
    geom_point(data = num.clust[1:i,],
                            color = 'red')+
    geom_line(data = num.clust[1:i,],
              color = 'red')
  
  
  p4 <- cowplot::plot_grid(p.a.i,
                           p.b.i,
                           p.c.i,ncol = 1,align = 'hv',axis = 'tblr') 
  p4
  
  # save plot as png
  setwd("/Users/msbr/Dropbox/a-s-omics/sam_sandbox")
  png(filename = paste(length(foi),'Features | Clusters | Entropy Summary Plot.png'),width = 18,height = 12,units = 'in',res = 300)
  print(cowplot::plot_grid(p1,p2,p3,p4,ncol = 2,align = 'hv',axis = 't'))
  dev.off()

  
  }

