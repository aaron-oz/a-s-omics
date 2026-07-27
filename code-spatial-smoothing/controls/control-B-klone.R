## Control B (R2.8) -- self-contained ckpt array task: ONE coordinate-shuffle
## permutation. Refits all genes' Poisson fields on spatially-permuted counts
## (paper's bru model), convolves the TRUE FANTOM pairings (common minimum),
## clusters (scale -> KNN/SNN/Louvain, no PCA), scores spatial coherence (ELSA),
## and appends a single summary line to $OUT in home. Idempotent + requeue-safe:
## overwrites its own perm row. Transients stay in node-local $TMPDIR.
##
## env: ASOMICS_ROOT (repo root on klone), OUT (summary csv in home),
##      SLURM_ARRAY_TASK_ID (permutation id), SLURM_CPUS_PER_TASK (parallel fits)
## perm_id 0 is the REAL (unpermuted) reference run.

suppressMessages({
  library(data.table); library(INLA); library(fmesher); library(splancs); library(inlabru)
  library(parallel); library(RcppHNSW); library(igraph); library(Matrix); library(elsa); library(raster)
})
root <- Sys.getenv("ASOMICS_ROOT", "."); setwd(root)
ctrl <- "data-outputs/mouse-embryo-raw/controls"
outcsv <- Sys.getenv("OUT", file.path(ctrl, "control-B-summary.csv"))
perm_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
ncore  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
NSAMP  <- as.integer(Sys.getenv("NSAMP", "1000"))   # posterior draws for the median field
inla.setOption(num.threads = "1:1")            # 1 thread/fit; parallelism via mclapply
tmpd <- Sys.getenv("TMPDIR", tempdir())

tic <- function(){ now <- proc.time(); function() round((proc.time()-now)[3],1) }

## ---- inputs: compact sub-region data, fixed mesh, FANTOM mechanisms ----
sub <- readRDS(file.path(ctrl, "bin50-subregion-raw.rds"))
locs <- unique(sub[, .(x,y)]); setkey(locs, x, y); nloc <- nrow(locs)
genes <- sort(unique(sub$feat)); GENE <- toupper(genes)
countwide <- dcast(sub, x + y ~ feat, value.var="feat.count"); setkey(countwide, x, y)
countwide <- countwide[locs]
numi <- unique(sub[, .(x,y,total.count)]); setkey(numi, x, y); numi <- numi[locs]$total.count
coord <- as.matrix(locs); rm(sub); gc()

load("data-inputs/mouse-embryo-raw/fantom-hierarchy.Robj")   # loads object 'output'
hier <- as.data.table(output)
mech <- unique(hier[, .(L=LIGAND,R=RECEPTOR)])[L %in% GENE & R %in% GENE]
mech[, li:=match(L,GENE)][, ri:=match(R,GENE)]

sb <- as.points(matrix(c(15000,25000,25000,15000,15000, 11500,11500,18000,18000,11500), ncol=2))
maxedge <- diff(range(sb[,1]))/(3*5); boundouter <- diff(range(sb[,1]))/3
mesh <- fm_mesh_2d(boundary=list(fm_segm(c(fm_segm(sb)))), max.edge=c(.15,2)*maxedge,
                   offset=c(maxedge, boundouter/1.5), cutoff=maxedge/10, min.angle=21)
matern.pri <- c(500,.95,.1,.05)
matern <- inla.spde2.pcmatern(mesh=mesh, alpha=2, constr=TRUE,
                              prior.range=matern.pri[1:2], prior.sigma=matern.pri[3:4])
bru_options_set(control.compute=list(cpo=FALSE))

fit_field <- function(cnt, tot){
  rd <- data.table(x=coord[,1], y=coord[,2], feat.count=cnt, total.count=tot)
  comps <- ~ fint(1) + ffield(cbind(x,y), model=matern)
  lik <- bru_obs(family="poisson", data=rd, formula=feat.count ~ fint+ffield, E=total.count)
  fit <- try(bru(comps, lik, options=list(bru_verbose=0, bru_max_iter=1)), silent=TRUE)
  if(inherits(fit,"try-error")) return(rep(NA_real_, nrow(rd)))
  pr <- try(predict(fit, rd, ~ exp(fint+ffield)*total.count, n.samples=NSAMP), silent=TRUE)
  if(inherits(pr,"try-error")) return(rep(NA_real_, nrow(rd)))
  pr$median
}

## ---- apply permutation (0 = real reference) ----
if(perm_id == 0){ perm <- seq_len(nloc) } else { set.seed(1000+perm_id); perm <- sample(nloc) }
numi.p <- numi[perm]

## ---- refit every gene in parallel ----
Tt <- tic()
cols <- lapply(genes, function(g) countwide[[g]][perm])
G <- do.call(cbind, mclapply(cols, function(cnt) fit_field(cnt, numi.p), mc.cores=ncore))
colnames(G) <- genes
message(sprintf("perm %d: refit %d genes on %d cores in %ss", perm_id, length(genes), ncore, Tt()))

## ---- convolve (true pairing) -> cluster -> ELSA ----
M <- pmin(G[,mech$li,drop=FALSE], G[,mech$ri,drop=FALSE])
Xs <- scale(M); Xs[is.na(Xs)] <- 0; Xs[Xs>10]<-10; Xs[Xs< -10]<- -10
idx <- hnsw_knn(Xs, k=20, distance="euclidean", n_threads=ncore)$idx
A <- sparseMatrix(i=rep(seq_len(nloc),each=20), j=as.vector(t(idx)), x=1, dims=c(nloc,nloc))
ov <- as(A %*% t(A), "TsparseMatrix"); jac <- ov@x/(40-ov@x); keep <- jac>(1/15) & ov@i!=ov@j
g <- simplify(graph_from_data_frame(data.frame(from=ov@i[keep]+1,to=ov@j[keep]+1,weight=jac[keep]),directed=FALSE))
memb <- rep(NA_integer_,nloc); memb[as.integer(V(g)$name)] <- membership(cluster_louvain(g, resolution=1))
if(any(is.na(memb))) memb[is.na(memb)] <- max(memb,na.rm=TRUE)+seq_len(sum(is.na(memb)))
em <- function(dd){ r <- rasterFromXYZ(cbind(coord[,1],coord[,2],as.integer(memb))); mean(elsa(r,d=dd,categorical=TRUE)@data@values,na.rm=TRUE) }
row <- data.table(perm_id=perm_id, n_clusters=length(unique(memb)), elsa_d50=em(50), elsa_d100=em(100))

## ---- one file per permutation (race-free for parallel array tasks) ----
## collate afterwards with: rbindlist(lapply(list.files(dir,'control-B-perm.*csv',full=T),fread))
outdir <- dirname(outcsv); dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
permfile <- file.path(outdir, sprintf("control-B-perm%03d.csv", perm_id))
fwrite(row, permfile)
message(sprintf("perm %d: clusters=%d ELSA d50=%.3f d100=%.3f -> %s",
                perm_id, row$n_clusters, row$elsa_d50, row$elsa_d100, permfile))
