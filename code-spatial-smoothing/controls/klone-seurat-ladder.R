.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
## Reproduce the paper's EXACT Seurat clustering (Sup7 recipe) on klone, then run the
## representation ladder on the same footing. Uses annoy (default, works on rocker;
## crashes only in the local distrobox). Validates rung4 == published 20 domains.
suppressMessages({library(Seurat); library(data.table); library(elsa); library(raster)})
root <- Sys.getenv("ASOMICS_ROOT","."); setwd(root); ctrl <- "stage"
GG<-readRDS(file.path(ctrl,"gene-field-matrix.rds")); G<-GG$G; xy<-GG$xy; GENE<-toupper(colnames(G))
load(file.path(ctrl,"fantom-hierarchy.Robj")); hier<-as.data.table(output)
mech<-unique(hier[,.(L=LIGAND,R=RECEPTOR)])[L%in%GENE & R%in%GENE]; mech[,li:=match(L,GENE)][,ri:=match(R,GENE)]
pub<-as.data.table(readRDS(file.path(ctrl,"paper-domain-labels.rds")))
ann<-fread(file.path(ctrl,"e16-fov-annotation.csv")); mp<-fread(file.path(ctrl,"mosta-plus-annotation.csv"))
key<-data.table(x=xy[,1],y=xy[,2]); pubv<-pub[key,on=c("x","y")]$domain; annv<-ann[key,on=c("x","y")]$annotation; mpv<-mp[key,on=c("x","y")]$ann_plus

## Sup7-exact: CreateAssayObject(counts) -> ScaleData(all features) -> FindNeighbors(features, no PCA) -> FindClusters res 1
seurat_cluster <- function(M, seed=42){
  Mt<-t(M); colnames(Mt)<-paste0("c",seq_len(ncol(Mt))); rownames(Mt)<-paste0("f",seq_len(nrow(Mt)))
  ass<-CreateAssayObject(counts=Mt); o<-CreateSeuratObject(ass, assay="RNA")
  o<-ScaleData(o, features=rownames(o), verbose=FALSE)
  o<-FindNeighbors(o, features=rownames(o), dims=NULL, verbose=FALSE)   # annoy (default)
  o<-FindClusters(o, resolution=1, random.seed=seed, verbose=FALSE)
  as.integer(o$seurat_clusters) }
em<-function(l){r<-rasterFromXYZ(cbind(xy[,1],xy[,2],as.integer(l)));mean(elsa(r,d=50,categorical=TRUE)@data@values,na.rm=TRUE)}
ARI<-function(a,b){t<-table(a,b);n<-sum(t);ci<-function(x)sum(choose(x,2));(ci(t)-ci(rowSums(t))*ci(colSums(t))/choose(n,2))/((ci(rowSums(t))+ci(colSums(t)))/2-ci(rowSums(t))*ci(colSums(t))/choose(n,2))}
pur<-function(l,g){t<-table(l,g);sum(apply(t,1,max))/sum(t)}

uni<-sort(unique(c(mech$li,mech$ri)))
reps<-list(rung2_union=G[,uni,drop=FALSE], rung3_product=G[,mech$li,drop=FALSE]*G[,mech$ri,drop=FALSE],
           rung4_commonmin=pmin(G[,mech$li,drop=FALSE],G[,mech$ri,drop=FALSE]))
out<-list()
for(nm in names(reps)){ t0<-proc.time()[3]; lab<-seurat_cluster(reps[[nm]])
  out[[nm]]<-data.table(rep=nm, n_feat=ncol(reps[[nm]]), n_clust=length(unique(lab)), elsa=em(lab),
    ARI_vs_published=ARI(lab,pubv), purity_anat=pur(lab,annv), purity_plus=pur(lab,mpv),
    homo_anat=NA_real_)
  saveRDS(lab, file.path(root,paste0("labels-",nm,".rds")))
  message(sprintf("%s: %d clusters in %.0fs; ARI_vs_published=%.3f", nm, length(unique(lab)), proc.time()[3]-t0, ARI(lab,pubv))) }
res<-rbindlist(out, fill=TRUE); print(res); fwrite(res, file.path(root,"klone-seurat-ladder.csv"))
cat("\nVALIDATION: rung4 ARI vs published should be ~1.0 if the Seurat recipe reproduces the paper.\n")
