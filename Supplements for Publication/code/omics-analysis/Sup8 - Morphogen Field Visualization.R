## 2025-07-30
## MSBR
## experimentation for figure 6

# globals
options(future.globals.maxSize = 10.0 * 1e9)


# wd
setwd("~/Dropbox/a-s-omics/sam_sandbox/morphogen fields")

#pack
require(elsa)
require(ggplot2)
require(tmap)
require(tidyr)
require(Seurat)

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

# load data
load('/Users/msbr/Dropbox/a-s-omics/sam_sandbox/clustering vs features/m.est.min/clustering vs features experiment/object.feature.depth.experiment.2025-07-12.Robj')
obj

# isolate assay of interest
meta <- obj@meta.data
ggplot(data = meta,
       aes(x=x,
           y=y,
           color = `1481featureClusters`))+
  geom_point(size = 0.5)+
  theme_classic()+DarkTheme()+scale_color_manual(values = cols.global)+ggtitle(paste('m.est.min','| Spatial'))

# identify features of interest
DefaultAssay(obj) <- 'm.est.min'
genes.to.study <- FindAllMarkers(obj,logfc.threshold = 1)

# identify range of interest
x.range <- c(18000,20000)
y.range <- c(16000,18000)
subset.meta <- meta[meta$x >= x.range[1] & meta$x <= x.range[2] & meta$y >= y.range[1] & meta$y <= y.range[2],]
subset.obj <- subset(obj,cells = rownames(subset.meta))
subset.data <- subset.obj@assays$m.est.min@layers$counts
rownames(subset.data) <- rownames(subset.obj)
colnames(subset.data) <- colnames(subset.obj)

# make plot showing ROI
ggplot(data = meta,
       aes(x=x,
           y=y,
           color = `1481featureClusters`))+
  geom_point(size = 3,shape = 15)+
  annotate(geom = "rect", 
           ymax = y.range[2], 
           ymin = y.range[1], 
           xmax = x.range[2], 
           xmin = x.range[1], 
           colour = "red",
           linetype = 'dashed',
           fill = NA)+
  theme_classic()+DarkTheme()+scale_color_manual(values = cols.global)+ggtitle(paste('m.est.min','| Spatial'))

# make plots showing FOI
goi <- c('SHH-PTCH1','RSPO1-LGR5','EFNA5-EPHA1','BMP8B-BMPR1A','LAMA5-ITGB1','SLIT1-ROBO1','IGF2-IGF1R','DCN-EGFR','FGF7-FGFR2','THBS2-ITGB1','DLK1-NOTCH2')
goi <- rownames(subset.obj)

to.plot <- cbind(subset.meta,t(subset.data[goi,]))
ggplot(to.plot,
       aes(x = x,
           y = y,
           color = `A2M-LRP1`))+
geom_point(shape = 15,size=25)

# 3d
library(plotly)

for (i in 1:length(goi)){
  temp <- to.plot[,c('x','y',goi[i])]
  temp <- reshape2::dcast(temp,x~y)
  rownames(temp) <- temp[,1]
  temp <- as.matrix(temp[,-1])
  temp2 <- t(temp)
  
  fig <- plot_ly(z = ~temp2, 
                 type = "surface") %>%
  #layout(scene = list(zaxis = list(range = c(0,1))))
  add_surface() %>% layout(
    scene = list(
      aspectmode = "manual", # Set aspectmode to "manual" for custom ratios
      aspectratio = list(x = 1, y = 1, z = 0.075), # Define desired ratios (e.g., z half of x and y)
      camera = list(
        eye = list(x = 1, y = -1, z = 1), # Adjust x, y, z for desired eye position
        center = list(x = 0, y = 0, z = 0), # Point the camera towards the origin
        up = list(x = 0, y = 0, z = 1) # Define the 'up' direction (Z-axis in this case)
      )
    ))
  fig
  htmlwidgets::saveWidget(as_widget(fig), paste(goi[i],".html",sep=''))
  
}

