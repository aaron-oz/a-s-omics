## 2025-07-26
## This script recreates the works from May 2025, but makes the rows unique in each Seurat object, which isn't naturally he case because of how we have indexed the ligand only and receptor only assays

## MSBR
options(future.globals.maxSize = 10.0 * 1e9)
# Set wd
setwd("~/Dropbox/a-s-omics/sam_sandbox/HD Embeddings")

# colors and make unique
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

cols.global <- unique(c(cols.epi,cols.endo,cols.mes,cols.lymphoid,cols.myeloid,'white','yellow','#A88FAC','#D2D68D','#CBD4C2','red','blue','purple','magenta'))

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

# turn each layer of information into a seurat object -> make rows unique -> pca -> cluster -> embed -> store
layer.names <- c("l.d","r.d","l.est","r.est","l.est.d","r.est.d","m.data.min",
                 "m.data.prod","m.data.gm","m.data.k.01","m.data.k.5","m.data.k1","m.data.d.min",
                 "m.data.d.prod" ,"m.data.d.gm","m.data.d.k.01", "m.data.d.k.5",  "m.data.d.k1",
                 "m.est.min","m.est.prod","m.est.gm","m.est.k.01","m.est.k.5","m.est.k1",
                 "m.est.d.min","m.est.d.prod","m.est.d.gm","m.est.d.k.01","m.est.d.k.5","m.est.d.k1",'l','r')

assay.list <- list()
for(k in 1:length(layer.names)){
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
  
  ### new 2025-07-26 MAKE ROWS UNIQUE (remove duplicate ligand or receptor info rows which arise in AOZ data due to fundamental signaling mechanism structure)
  ligand.set <- c("l.d","l.est","l.est.d",'l')
  receptor.set <- c("r.d","r.est","r.est.d",'r')
  
  if(layer.names[k] %in% ligand.set){
    temp.2 <- as.matrix(temp)
    rownames(temp.2) <- stringr::str_split_fixed(rownames(temp.2),pattern='_',2)[,1]
    temp.2 <- temp.2[unique(rownames(temp.2)),]
    temp <- temp.2
  }
  if(layer.names[k] %in% receptor.set){
    temp.2 <- as.matrix(temp)
    rownames(temp.2) <- stringr::str_split_fixed(rownames(temp.2),pattern='_',2)[,2]
    temp.2 <- temp.2[unique(rownames(temp.2)),]
    temp <- temp.2
  }

  
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
  setwd("/Users/msbr/Dropbox/a-s-omics/sam_sandbox/HD Embeddings")
  png(filename = paste(layer.name,'HD Embeddings.png'),width = 18,height = 6,units = 'in',res = 300)
  print(p1|p2)
  dev.off()
  # name and save
  names(assay.list) <- layer.names[1:length(assay.list)]
  setwd("/Users/msbr/Dropbox/a-s-omics/sam_sandbox/HD Embeddings")
  #save(assay.list,file = 'assay.list.2025-07-26.Robj')
  
  save(temp,file = paste(layer.name,'assay.2025-07-25.Robj',sep=''))
  
  message(paste('Computation complete for Layer #',k,': ','"',layer.name,'"',sep=''))
}

save(assay.list,file = 'assay.list.2025-07-26.Robj')
