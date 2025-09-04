### ELSA work, updated 2025-07-27

## INCLUDES CONVOLUTION PLOTS AT THE END

## MSBR

# globals
options(future.globals.maxSize = 10.0 * 1e9)


# wd
setwd("~/Dropbox/a-s-omics/sam_sandbox/ELSA/ELSA 2nd run (July 2025)")

#pack
require(elsa)
require(ggplot2)
require(tmap)
require(tidyr)

## colors
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

# elsa https://cran.r-project.org/web/packages/elsa/vignettes/elsa.html
file <- system.file('external/lc_example.grd',package='elsa') 
lc <- raster(file)
plot(lc,main='Land cover: a categorical map')
elc <- elsa(lc, d=2000,categorical = T)
cl <- colorRampPalette(c('darkblue','yellow','red','black'))(100) # specifying a color scheme
plot(elc, col=cl, main='ELSA')

# make a raster object from each seurat dataset metadata...https://gis.stackexchange.com/questions/404443/create-raster-from-dataframe-containing-characters-as-values-in-r

load("~/Dropbox/a-s-omics/sam_sandbox/HD Embeddings/assay.list.2025-07-26.Robj")
names(assay.list)

assay.names.list <- names(assay.list)

assay.raw.elsa.values <- list()
assay.mean.elsa.trends <- list()
for(i in 1:length(assay.names.list)){
  print(i)
  assay.of.interest <- assay.names.list[i]
  
  #isolate assay of interest
  temp <- assay.list[[assay.of.interest]]@meta.data
  
  # temp$seurat_clusters <- as.numeric(temp$seurat_clusters)
  # temp$seurat_clusters <- temp$seurat_clusters+1
  # temp$seurat_clusters <- as.character(temp$seurat_clusters)
  temp$seurat_clusters <- paste('Cluster_',temp$seurat_clusters,sep='')
  temp$seurat_clusters <- as.factor(temp$seurat_clusters )
  temp$seurat_clusters_numeric <- as.numeric(temp$seurat_clusters)
  
  # convert to raster
  r <- rasterFromXYZ(temp[,c("x","y","seurat_clusters_numeric")])
  r[] <- factor(levels(temp$seurat_clusters)[r[]])
  
  # print a quick view of the cluster structure
  png(filename = paste(assay.of.interest,'Cluster tmap.png'),width = 7,height = 5,res = 300,units = 'in')
  print(tm_shape(r) + tm_raster())
  dev.off()
  
  # run elsa for a variety of distance parameters
  dist.list <-seq(50,1000,50)
  elsa.by.dist <- list()
  for(j in 1:length(dist.list)){
    test <- elsa(r, d=dist.list[j],categorical = T)
    
    png(filename = paste(assay.of.interest,' ELSA | d = ',dist.list[j],' | Enterogram', '.png',sep = ''),width = 7,height = 5,res = 300,units = 'in')
    print(plot(test, col=cl, main='ELSA'))
    dev.off()
    
    elsa.by.dist[[j]] <- test@data@values
    
    png(filename = paste(assay.of.interest,' ELSA | d = ',dist.list[j],' | Distribution', '.png',sep = ''),width = 7,height = 5,res = 300,units = 'in')
    print(ggplot(data = data.frame(elsa = elsa.by.dist[[j]]),
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
    mean.elsa[[k]] <- mean(elsa.by.dist[[k]])
  }
  # Stash mean elsa values (trend with increasing d)
  assay.mean.elsa.trends[[i]] <- unlist(mean.elsa)
}

# name outputs for each assay
names(assay.raw.elsa.values) <- assay.names.list[1:length(assay.raw.elsa.values)]
names(assay.mean.elsa.trends) <- assay.names.list[1:length(assay.mean.elsa.trends)]
names(assay.raw.elsa.values)
names(assay.mean.elsa.trends)

# save for later
save(assay.raw.elsa.values,file = 'assay.raw.elsa.values.2025-07-27.Robj')
save(assay.mean.elsa.trends,file = 'assay.mean.elsa.trends.2025-07-27.Robj')

# Plottting
## MSBR revisited 2025-07-27

to.plot <- data.frame(assay.mean.elsa.trends)
dist.list <-seq(50,1000,50)
to.plot$param_d <- dist.list
to.plot <- pivot_longer(to.plot,
                        cols = c('l.d':'r'),
                        names_to = 'assay',
                        values_to = 'elsa')
#View(to.plot)
assays.to.include <- c("l.est", "r.est",
                       "m.data.min","m.data.prod","m.data.gm","m.data.k.01","m.data.k.5","m.data.k1",
                       "m.est.min","m.est.prod","m.est.gm","m.est.k.01","m.est.k.5","m.est.k1",
                       "l","r")
assays.to.include.2 <- c("l.est", "r.est",
                       "m.est.min","m.est.prod","m.est.gm",
                       "l","r")

assay.cols <- c('red','blue',
                'yellow','green','darkgrey','white',"violetred",'darkgreen',
                'orange','purple',"lightseagreen",
                'lightgreen',"lightpink","navyblue",
                "darkred","steelblue")
names(assay.cols) <- assays.to.include

pdf(file = 'elsa.v.paramD.pdf',width = 7, height = 3.5)
ggplot(to.plot[to.plot$assay %in% assays.to.include.2,],
       aes(x = param_d,
           y = elsa,
           color = assay,
           group = assay))+
  geom_point()+
  geom_line()+
  # geom_label(data = to.plot[to.plot$param_d == max(to.plot$param_d),],
  #   aes(label = assay), nudge_x = 0.35, size = 4)+
  theme_classic()+
  ylim(0,NA)+
  Seurat::DarkTheme()+
  scale_color_manual(values = assay.cols)
dev.off()

# scatter against number of clusters?

num.clust <- list()
for(i in 1:length(assay.list)){
  num.clust[[i]] <- length(names(table(assay.list[[i]]@meta.data$seurat_clusters)))
}
names(num.clust) <- names(assay.list)
num.clust
to.plot[to.plot$param_d == 100,]

scatter.dat <- cbind(data.frame(num.clust = unlist(num.clust)),to.plot[to.plot$param_d == 100,])

scatter.dat$ratio_nClust.ELSA <- scatter.dat$num.clust/scatter.dat$elsa

save(scatter.dat,file = 'scatter.dat.2025-08-07.Robj')

# ggplot(scatter.dat[scatter.dat$assay %in% assays.to.include,],
#        aes(x = num.clust,
#            y = elsa,
#            color = assay))+
#   geom_point(size=4)+
#   Seurat::DarkTheme()+
#   scale_color_manual(values = cols.global[sample(length(cols.global))])
# 


# assay.cols <- c("darkred",'red',
#                 'purple','lightgreen','orange',
#                 "lightseagreen","navyblue",
#                 'blue','yellow','green','grey','white',
#                 "steelblue",  "violetred",'darkgreen',"lightpink") 

                #            "royalblue"       "plum4"           "lightgoldenrod"         "dimgray"         "deeppink"       
                #  [21] "red2"            "paleturquoise1"  "palevioletred"   "#E0CA3C"         "#989FCE"         "#261132"         "#F37748"         "#AFFC41"         "#347FC4"         "#0B5351"        
                #  [31] "#00A9A5"         "#4E8098"         "#D81159"         "#AB3428"         "#01BAEF"         "#F42C04"         "#610345"         "#044B7F"         "orchid4"         "purple4"        
                #  [41] "plum1"           "olivedrab2"          "mediumvioletred" "sienna"          "orange"              "lightseagreen"   "mediumpurple4"   "#F8F4A6"        
                # "yellow"          "#A88FAC"         "#D2D68D"         "#CBD4C2"         "red"             "blue"            "magenta"        )
                # 

pdf(file = 'elsa.v.numclust.pdf',width = 7, height = 3.5)
ggplot(scatter.dat[scatter.dat$assay %in% assays.to.include,],
       aes(x = num.clust,
           y = elsa,
           color = assay))+
  geom_point(size=4)+
  ylim(0,NA)+
  Seurat::DarkTheme()+
  scale_color_manual(values = assay.cols)
dev.off()

#### plots to communicate convolution comparisons [could move to a dedicated script?]
## MSBR 2025-07-28

ConvolutionPlot <- function(convolved.assay = 'm.est.min',
                            ligand.assay = 'l.est',
                            receptor.assay = 'r.est',
                            ligand = 'WNT4',
                            receptor = 'FZD6',
                            viridis.option = 'A',
                            color.thresh = 0.25,
                            make.pdf = T,
                            pdf.width = 9,
                            pdf.height = 3.75){
  # collect info to use
  convolved.info <- assay.list[[convolved.assay]]
  convolved.counts <- convolved.info@assays[[convolved.assay]]@layers$counts
  rownames(convolved.counts) <- rownames(convolved.info)
  colnames(convolved.counts) <- colnames(convolved.info)
  
  ligand.info <- assay.list[[ligand.assay]]
  ligand.counts <- ligand.info@assays[[ligand.assay]]@layers$counts
  rownames(ligand.counts) <- rownames(ligand.info)
  colnames(ligand.counts) <- colnames(ligand.info)
  
  receptor.info <- assay.list[[receptor.assay]]
  receptor.counts <- receptor.info@assays[[receptor.assay]]@layers$counts
  rownames(receptor.counts) <- rownames(receptor.info)
  colnames(receptor.counts) <- colnames(receptor.info)
  
  # make to.plots
  mech <- paste(ligand,receptor,sep='-')
  to.plot.conv <- cbind(convolved.info@meta.data,
                        mech = convolved.counts[mech,])
  to.plot.lig <- cbind(ligand.info@meta.data,
                        ligand = ligand.counts[ligand,])
  to.plot.rec <- cbind(ligand.info@meta.data,
                       receptor = receptor.counts[receptor,])
  
  # plot
p1 <-  ggplot(to.plot.lig,
         aes(x=x,
             y=y,
             color = ligand))+
    geom_point(shape=15)+
    scale_color_viridis_c(option=viridis.option,
                          limits = c(0,color.thresh*max(to.plot.lig$ligand)),
                          oob = scales::squish)+
  ggtitle(paste(ligand,'|',ligand.assay,sep=' '))+
  Seurat::DarkTheme()+
  Seurat::NoAxes()+theme(legend.title=element_blank())+  theme(axis.line = element_line(colour = "black"),
                                                               panel.grid.major = element_blank(),
                                                               panel.grid.minor = element_blank(),
                                                               panel.border = element_blank(),
                                                               panel.background = element_blank())

p2 <-  ggplot(to.plot.rec,
         aes(x=x,
             y=y,
             color = receptor))+
    geom_point(shape=15)+
    scale_color_viridis_c(option=viridis.option,
                          limits = c(0,color.thresh*max(to.plot.rec$receptor)),
                          oob = scales::squish)+
  ggtitle(paste(receptor,'|',receptor.assay,sep=' '))+
  Seurat::DarkTheme()+
  Seurat::NoAxes()+theme(legend.title=element_blank())+  theme(axis.line = element_line(colour = "black"),
                                                               panel.grid.major = element_blank(),
                                                               panel.grid.minor = element_blank(),
                                                               panel.border = element_blank(),
                                                               panel.background = element_blank())
  
p3 <-  ggplot(to.plot.conv,
         aes(x=x,
             y=y,
             color = mech))+
    geom_point(shape=15)+
    scale_color_viridis_c(option=viridis.option,
                          limits = c(0,color.thresh*max(to.plot.conv$mech)),
                          oob = scales::squish)+
  ggtitle(paste(mech,'|',convolved.assay,sep=' '))+
  Seurat::DarkTheme()+
  Seurat::NoAxes()+theme(legend.title=element_blank())+  theme(axis.line = element_line(colour = "black"),
                                                               panel.grid.major = element_blank(),
                                                               panel.grid.minor = element_blank(),
                                                               panel.border = element_blank(),
                                                               panel.background = element_blank())

p4 <- cowplot::plot_grid(p1,p2,ncol = 1,align = 'hv')
cowplot::plot_grid(p4,p3,ncol = 2,rel_widths = c(1,2))
plot.out <- cowplot::plot_grid(p4,p3,ncol = 2,rel_widths = c(1,2))
if(make.pdf == T){
  pdf(file = paste(mech,'|',convolved.assay,'.pdf',sep=' '),width = pdf.width,height=pdf.height)
  print(plot.out)
  dev.off()
}
return(plot.out)
}


ConvolutionPlot(convolved.assay = 'm.est.min',
                ligand.assay = 'l.est',
                receptor.assay = 'r.est',
                ligand = 'WNT4',
                receptor = 'FZD6',
                viridis.option = 'A',
                color.thresh = 0.2)
ConvolutionPlot(convolved.assay = 'm.est.gm',
                ligand.assay = 'l.est',
                receptor.assay = 'r.est',
                ligand = 'WNT4',
                receptor = 'FZD6',
                viridis.option = 'A',
                color.thresh = 0.2)
ConvolutionPlot(convolved.assay = 'm.est.prod',
                ligand.assay = 'l.est',
                receptor.assay = 'r.est',
                ligand = 'WNT4',
                receptor = 'FZD6',
                viridis.option = 'A',
                color.thresh = 0.2)

ConvolutionPlot(convolved.assay = 'm.est.min',
                ligand.assay = 'l.est',
                receptor.assay = 'r.est',
                ligand = 'RSPO1',
                receptor = 'LGR6',
                viridis.option = 'A',
                color.thresh = 0.5)
ConvolutionPlot(convolved.assay = 'm.est.gm',
                ligand.assay = 'l.est',
                receptor.assay = 'r.est',
                ligand = 'RSPO1',
                receptor = 'LGR6',
                viridis.option = 'A',
                color.thresh = 0.5)
ConvolutionPlot(convolved.assay = 'm.est.prod',
                ligand.assay = 'l.est',
                receptor.assay = 'r.est',
                ligand = 'RSPO1',
                receptor = 'LGR6',
                viridis.option = 'A',
                color.thresh = 0.5)

ConvolutionPlot(convolved.assay = 'm.est.min',
                ligand.assay = 'l.est',
                receptor.assay = 'r.est',
                ligand = 'DHH',
                receptor = 'PTCH1',
                viridis.option = 'A',
                color.thresh = 0.5)
ConvolutionPlot(convolved.assay = 'm.est.gm',
                ligand.assay = 'l.est',
                receptor.assay = 'r.est',
                ligand = 'DHH',
                receptor = 'PTCH1',
                viridis.option = 'A',
                color.thresh = 0.5)
ConvolutionPlot(convolved.assay = 'm.est.prod',
                ligand.assay = 'l.est',
                receptor.assay = 'r.est',
                ligand = 'DHH',
                receptor = 'PTCH1',
                viridis.option = 'A',
                color.thresh = 0.5)

ConvolutionPlot(convolved.assay = 'm.est.min',
                ligand.assay = 'l.est',
                receptor.assay = 'r.est',
                ligand = 'VEGFA',
                receptor = 'KDR',
                viridis.option = 'A',
                color.thresh = 0.5)
ConvolutionPlot(convolved.assay = 'm.est.gm',
                ligand.assay = 'l.est',
                receptor.assay = 'r.est',
                ligand = 'VEGFA',
                receptor = 'KDR',
                viridis.option = 'A',
                color.thresh = 0.5)
ConvolutionPlot(convolved.assay = 'm.est.prod',
                ligand.assay = 'l.est',
                receptor.assay = 'r.est',
                ligand = 'VEGFA',
                receptor = 'KDR',
                viridis.option = 'A',
                color.thresh = 0.5)



