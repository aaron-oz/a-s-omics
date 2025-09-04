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
# set up
require(Seurat)
require(ggplot2)
file.path.input <- "/Users/msbr/Dropbox/a-s-omics/data-outputs/mouse-embryo-raw/2025-04-10/convolution-data"
file.path.list <- list.files(file.path.input, full.names = T)
file.list <- list.files(file.path.input)

# load in data from aoz [takes a few minutes total]
data.list <- list()
for(i in 1:length(file.list)){
  print(i)
  data.list[[i]] <- read.csv(file.path.list[i])
}

# name each mechanism
mech.names <- gsub(".csv", "", file.list)
mech.names <- gsub("convolution-calcs-", "", mech.names)
names(data.list) <- mech.names

# turn each layer of information into a seurat object -> pca -> cluster -> embed -> store
layer.names <- c("l.d","r.d","l.est","r.est","l.est.d","r.est.d","m.data.min",
                 "m.data.prod","m.data.gm","m.data.k.01","m.data.k.5","m.data.k1","m.data.d.min",
                 "m.data.d.prod" ,"m.data.d.gm","m.data.d.k.01", "m.data.d.k.5",  "m.data.d.k1",
                 "m.est.min","m.est.prod","m.est.gm","m.est.k.01","m.est.k.5","m.est.k1",
                 "m.est.d.min","m.est.d.prod","m.est.d.gm","m.est.d.k.01","m.est.d.k.5","m.est.d.k1",'l','r')
assay.list <- list()
for(k in 27:length(layer.names)){
  gc()
  print(k)
  # define layer name to pull
  layer.name <- layer.names[k]
  print(layer.name)
  # isolate one assay
  temp <- data.frame()
  feature.list <- c()
  for(i in 1:length(data.list)){
    print(i)
    if(sum(data.list[[i]][[layer.name]])>0){
    temp <- rbind(temp,data.list[[i]][[layer.name]])
    feature.list <- c(feature.list,names(data.list)[i])
    }
  }
  #rownames(temp) <- names(data.list)
  rownames(temp) <- feature.list
  colnames(temp) <- paste(data.list[[1]]$x,data.list[[1]]$y,sep='.')
  #View(temp)
  # turn into seurat
  temp <- CreateAssay5Object(counts = Matrix::Matrix(as.matrix(temp),sparse = T),
                             data = Matrix::Matrix(as.matrix(temp),sparse = T))
  temp <- CreateSeuratObject(temp,
                             assay = layer.name)
  temp@meta.data$x <- data.list[[1]]$x
  temp@meta.data$y <- data.list[[1]]$y
  temp
  # scale
  temp <- ScaleData(temp)
  # variable features
  #temp <- FindVariableFeatures(temp,selection.method = 'disp')
  temp <- FindVariableFeatures(temp)
  # PCA
  temp <- RunPCA(temp,npcs = 100)
  ElbowPlot(temp,ndims = 100)
  DimHeatmap(temp,dims = 1:9,balanced = T,cells = 300)
  DimHeatmap(temp,dims = 10:18,balanced = T,cells = 300)
  DimHeatmap(temp,dims = 19:27,balanced = T,cells = 300)
  DimHeatmap(temp,dims = 28:36,balanced = T,cells = 300)
  DimHeatmap(temp,dims = 37:45,balanced = T,cells = 300)
  DimHeatmap(temp,dims = 46:54,balanced = T,cells = 300)
  DimHeatmap(temp,dims = 55:63,balanced = T,cells = 300)
  DimHeatmap(temp,dims = 64:72,balanced = T,cells = 300)
  DimHeatmap(temp,dims = 73:81,balanced = T,cells = 300)
  DimHeatmap(temp,dims = 82:90,balanced = T,cells = 300)
  DimHeatmap(temp,dims = 91:99,balanced = T,cells = 300)
  
  temp <- RunUMAP(temp,dims = 1:100)
  temp <- FindNeighbors(temp,dims = 1:100)
  temp <- FindClusters(temp,resolution = 1)
  
  # stash
  assay.list[[k]] <- temp
  
  # plot for fun
  to.plot <- data.frame(cluster = temp$seurat_clusters,
                        x = temp$x,
                        y = temp$y)
  p1 <- ggplot(to.plot,
               aes(x=x,
                   y=y,
                   colour = cluster))+
    geom_point(size = 0.5)+
    theme_classic()+DarkTheme()+scale_color_manual(values = cols.global)+ggtitle(paste(layer.name,'| Spatial'))
  p2 <- DimPlot(temp,label=T,label.color = 'white',cols = cols.global)+theme_classic()+DarkTheme()+NoLegend()+ggtitle(paste(layer.name,'| UMAP'))
  print(p1|p2)
  setwd("/Users/msbr/Dropbox/a-s-omics/sam_sandbox")
  png(filename = paste(layer.name,'HD Embeddings.png'),width = 18,height = 6,units = 'in',res = 300)
  print(p1|p2)
  dev.off()
  # name and save
  names(assay.list) <- layer.names[1:length(assay.list)]
  setwd("/Users/msbr/Dropbox/a-s-omics/sam_sandbox")
  #save(assay.list,file = 'assay.list.2025-05-16.Robj')
  
  message(paste('Computation complete for Layer #',k,': ','"',layer.name,'"',sep=''))
 }

save(assay.list,file = 'assay.list.2025-05-16.Robj')

# # name and save
# names(assay.list) <- layer.names[1:length(assay.list)]
# setwd("/Users/msbr/Dropbox/a-s-omics/sam_sandbox")
# save(assay.list,file = 'assay.list.2025-05-16.Robj')










####### legacy below #####

# isolate one assay
m.est.min <- data.frame()

for(i in 1:length(data.list)){
  print(i)
  m.est.min <- rbind(m.est.min,data.list[[i]]$m.est.min)
}
rownames(m.est.min) <- names(data.list)
colnames(m.est.min) <- paste(data.list[[1]]$x,data.list[[1]]$y,sep='.')
View(m.est.min)

# turn into seurat
test <- CreateAssay5Object(counts = m.est.min,
                           data = m.est.min)
test <- CreateSeuratObject(test,
                           assay = 'm.est.min')
test@meta.data$x <- data.list[[1]]$x
test@meta.data$y <- data.list[[1]]$y
test
# scale
test <- ScaleData(test)
# variable features
test <- FindVariableFeatures(test,selection.method = 'disp')
# PCA
test <- RunPCA(test,npcs = 100)
ElbowPlot(test,ndims = 100)
DimHeatmap(test,dims = 1:9,balanced = T,cells = 300)
DimHeatmap(test,dims = 10:18,balanced = T,cells = 300)
DimHeatmap(test,dims = 19:27,balanced = T,cells = 300)
DimHeatmap(test,dims = 28:36,balanced = T,cells = 300)
DimHeatmap(test,dims = 37:45,balanced = T,cells = 300)
DimHeatmap(test,dims = 46:54,balanced = T,cells = 300)
DimHeatmap(test,dims = 55:63,balanced = T,cells = 300)
DimHeatmap(test,dims = 64:72,balanced = T,cells = 300)
DimHeatmap(test,dims = 73:81,balanced = T,cells = 300)
DimHeatmap(test,dims = 82:90,balanced = T,cells = 300)
DimHeatmap(test,dims = 91:99,balanced = T,cells = 300)

test <- RunUMAP(test,dims = 1:100)

test <- FindNeighbors(test,dims = 1:100)

test <- FindClusters(test,resolution = 1)

# plot
to.plot <- data.frame(cluster = test$seurat_clusters,
                      x = test$x,
                      y = test$y)
p1 <- ggplot(to.plot,
             aes(x=x,
                 y=y,
                 colour = cluster))+
  geom_point(size = 0.5)+
  theme_classic()+DarkTheme()+scale_color_manual(values = cols.global)+ggtitle('m.est.min | Spatial')
p2 <- DimPlot(test,label=T,label.color = 'white',cols = cols.global)+theme_classic()+DarkTheme()+NoLegend()+ggtitle('m.est.min | UMAP')
p1|p2


# isolate one assay
m.est.k1 <- data.frame()

for(i in 1:length(data.list)){
  print(i)
  m.est.k1 <- rbind(m.est.k1,data.list[[i]]$m.est.k1)
}
rownames(m.est.k1) <- names(data.list)
colnames(m.est.k1) <- paste(data.list[[1]]$x,data.list[[1]]$y,sep='.')
View(m.est.k1)

# turn into seurat
test.2 <- CreateAssay5Object(counts = m.est.k1,
                           data = m.est.k1)
test.2 <- CreateSeuratObject(test.2,
                           assay = 'm.est.k1')
test.2@meta.data$x <- data.list[[1]]$x
test.2@meta.data$y <- data.list[[1]]$y
test.2
# scale
test.2 <- ScaleData(test.2)
# variable features
test.2 <- FindVariableFeatures(test.2,selection.method = 'disp')
# PCA
test.2 <- RunPCA(test.2,npcs = 100)
ElbowPlot(test.2,ndims = 100)
DimHeatmap(test.2,dims = 1:9,balanced = T,cells = 300)
DimHeatmap(test.2,dims = 10:18,balanced = T,cells = 300)
DimHeatmap(test.2,dims = 19:27,balanced = T,cells = 300)
DimHeatmap(test.2,dims = 28:36,balanced = T,cells = 300)
DimHeatmap(test.2,dims = 37:45,balanced = T,cells = 300)
DimHeatmap(test.2,dims = 46:54,balanced = T,cells = 300)
DimHeatmap(test.2,dims = 55:63,balanced = T,cells = 300)
DimHeatmap(test.2,dims = 64:72,balanced = T,cells = 300)
DimHeatmap(test.2,dims = 73:81,balanced = T,cells = 300)
DimHeatmap(test.2,dims = 82:90,balanced = T,cells = 300)
DimHeatmap(test.2,dims = 91:99,balanced = T,cells = 300)

test.2 <- RunUMAP(test.2,dims = 1:100)

test.2 <- FindNeighbors(test.2,dims = 1:100)

test.2 <- FindClusters(test.2,resolution = 1)

# plot
to.plot <- data.frame(cluster = test.2$seurat_clusters,
                      x = test.2$x,
                      y = test.2$y)
p3 <- ggplot(to.plot,
             aes(x=x,
                 y=y,
                 colour = cluster))+
  geom_point(size = 0.5)+
  theme_classic()+DarkTheme()+scale_color_manual(values = cols.global)+ggtitle('m.est.k1 | Spatial')
p4 <- DimPlot(test.2,label=T,label.color = 'white',cols = cols.global)+theme_classic()+DarkTheme()+NoLegend()+ggtitle('m.est.k1 | UMAP')
p3|p4
cowplot::plot_grid(p1,p2,p3,p4,nrow=2)


# isolate one assay
m.est.gm <- data.frame()

for(i in 1:length(data.list)){
  print(i)
  m.est.gm <- rbind(m.est.gm,data.list[[i]]$m.est.gm)
}
rownames(m.est.gm) <- names(data.list)
colnames(m.est.gm) <- paste(data.list[[1]]$x,data.list[[1]]$y,sep='.')
View(m.est.gm)

# turn into seurat
test.3 <- CreateAssay5Object(counts = m.est.gm,
                             data = m.est.gm)
test.3 <- CreateSeuratObject(test.3,
                             assay = 'm.est.gm')
test.3@meta.data$x <- data.list[[1]]$x
test.3@meta.data$y <- data.list[[1]]$y
test.3
# scale
test.3 <- ScaleData(test.3)
# variable features
test.3 <- FindVariableFeatures(test.3,selection.method = 'disp')
# PCA
test.3 <- RunPCA(test.3,npcs = 100)
ElbowPlot(test.3,ndims = 100)
DimHeatmap(test.3,dims = 1:9,balanced = T,cells = 300)
DimHeatmap(test.3,dims = 10:18,balanced = T,cells = 300)
DimHeatmap(test.3,dims = 19:27,balanced = T,cells = 300)
DimHeatmap(test.3,dims = 28:36,balanced = T,cells = 300)
DimHeatmap(test.3,dims = 37:45,balanced = T,cells = 300)
DimHeatmap(test.3,dims = 46:54,balanced = T,cells = 300)
DimHeatmap(test.3,dims = 55:63,balanced = T,cells = 300)
DimHeatmap(test.3,dims = 64:72,balanced = T,cells = 300)
DimHeatmap(test.3,dims = 73:81,balanced = T,cells = 300)
DimHeatmap(test.3,dims = 82:90,balanced = T,cells = 300)
DimHeatmap(test.3,dims = 91:99,balanced = T,cells = 300)

test.3 <- RunUMAP(test.3,dims = 1:100)

test.3 <- FindNeighbors(test.3,dims = 1:100)

test.3 <- FindClusters(test.3,resolution = 1)

# plot
to.plot <- data.frame(cluster = test.3$seurat_clusters,
                      x = test.3$x,
                      y = test.3$y)
p5 <- ggplot(to.plot,
             aes(x=x,
                 y=y,
                 colour = cluster))+
  geom_point(size = 0.5)+
  theme_classic()+DarkTheme()+scale_color_manual(values = cols.global)+ggtitle('m.est.gm | Spatial')
p6 <- DimPlot(test.3,label=T,label.color = 'white',cols = cols.global)+theme_classic()+DarkTheme()+NoLegend()+ggtitle('m.est.gm | UMAP')
p5|p6
cowplot::plot_grid(p1,p2,p3,p4,p5,p6,nrow=3)


# isolate one assay
m.data.min <- data.frame()

for(i in 1:length(data.list)){
  print(i)
  m.data.min <- rbind(m.data.min,data.list[[i]]$m.data.min)
}
rownames(m.data.min) <- names(data.list)
colnames(m.data.min) <- paste(data.list[[1]]$x,data.list[[1]]$y,sep='.')
View(m.data.min)

# turn into seurat
test.4 <- CreateAssay5Object(counts = m.data.min,
                             data = m.data.min)
test.4 <- CreateSeuratObject(test.4,
                             assay = 'm.data.min')
test.4@meta.data$x <- data.list[[1]]$x
test.4@meta.data$y <- data.list[[1]]$y
test.4
# scale
test.4 <- ScaleData(test.4)
# variable features
test.4 <- FindVariableFeatures(test.4,selection.method = 'disp')
# PCA
test.4 <- RunPCA(test.4,npcs = 100)
ElbowPlot(test.4,ndims = 100)
DimHeatmap(test.4,dims = 1:9,balanced = T,cells = 300)
DimHeatmap(test.4,dims = 10:18,balanced = T,cells = 300)
DimHeatmap(test.4,dims = 19:27,balanced = T,cells = 300)
DimHeatmap(test.4,dims = 28:36,balanced = T,cells = 300)
DimHeatmap(test.4,dims = 37:45,balanced = T,cells = 300)
DimHeatmap(test.4,dims = 46:54,balanced = T,cells = 300)
DimHeatmap(test.4,dims = 55:63,balanced = T,cells = 300)
DimHeatmap(test.4,dims = 64:72,balanced = T,cells = 300)
DimHeatmap(test.4,dims = 73:81,balanced = T,cells = 300)
DimHeatmap(test.4,dims = 82:90,balanced = T,cells = 300)
DimHeatmap(test.4,dims = 91:99,balanced = T,cells = 300)

test.4 <- RunUMAP(test.4,dims = 1:100)

test.4 <- FindNeighbors(test.4,dims = 1:100)

test.4 <- FindClusters(test.4,resolution = 1)

# plot
to.plot <- data.frame(cluster = test.4$seurat_clusters,
                      x = test.4$x,
                      y = test.4$y)
p7 <- ggplot(to.plot,
             aes(x=x,
                 y=y,
                 colour = cluster))+
  geom_point(size = 0.5)+
  theme_classic()+DarkTheme()+scale_color_manual(values = cols.global)+ggtitle('m.data.min | Spatial')
p8 <- DimPlot(test.4,label=T,label.color = 'white',cols = cols.global)+theme_classic()+DarkTheme()+NoLegend()+ggtitle('m.data.min | UMAP')
p7|p8
cowplot::plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,nrow=4)

# save for later
m.est.min.seur <- test
m.est.k1.seur <- test.2
m.est.gm.seur <- test.3
m.data.min.seur <- test.4
to.save <- list(test,
     test.2,
     test.3,
     test.4)
setwd("/Users/msbr/Dropbox/a-s-omics/sam_sandbox")
save(to.save,file = 'to.save.2025-05-15.Robj')

# experimental plotting
rownames(test@assays$m.est.min@layers$counts) <- rownames(test)
colnames(test@assays$m.est.min@layers$counts) <- colnames(test)
range(test$x)
range(test$y)
x.coord <- 15025
y.coord <- 17975
object <- test
plot.sig <- function(x.coord,
                     y.coord,
                     object){
  location <- paste(x.coord,y.coord,sep='.')
  to.plot <- object@assays$m.est.min@layers$counts[,location]
  to.plot <- data.frame(feature = rownames(object),
                        value = to.plot)
  #View(to.plot)
  ggplot(to.plot,
         aes(x = feature,
             y = value))+
    geom_bar(stat='identity')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    ggtitle(paste('Ligand-Receptor Signature for Location',location))+
    ylab('Connectivity')
  
}
range(test$x)
range(test$y)
x.coord <- 15025
y.coord <- 17675
plot.sig(x.coord = x.coord,
         y.coord = y.coord,
         object = object)
