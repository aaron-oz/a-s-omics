## this script launches  our model for different features. it can easily be parallelized


##########################################################
## set user selected params, in/out dirs, and load pkgs ##
##########################################################
repo.fp <- "path/to/repo"
repo.fp <- "~/Dropbox/genetics/a-s-omics/"
setwd(repo.fp)
source("./Supplements for Publication/code/spatial-modeling/supp-3-raw-v-norm-joint-feat-total-model/00-load-pkg-set-io-set-params.R")

poi.mod <- TRUE
zip.mod <- TRUE
zap.mod <- TRUE

if(poi.mod & zip.mod & zap.mod){
  run.lik.comp <- TRUE
}else{
  run.lik.comp <- FALSE
}

#######################
## load data, subset ##
#######################
## agg size set in set-params
load(file.path(i.d, glue("pre-fit-obj-binned-{agg.size}.Rdata")))
## subset, make boundary domain
if(!is.null(sub.bound)){
  in.sb <- which(splancs::inout(all.dat[, .(x, y)], sub.bound))
  sub.dat <- all.dat[in.sb, ]
  boundary.coords <- sub.bound
}else{
  sub.dat <- all.dat
  boundary.coords <- NULL #TODO
}
domain <- boundary.coords |> fmesher::fm_segm()

## subset to hierarchy mechanisms
sub.dat <- sub.dat[feat %in% f.to.mod, ]

########################################
## make mesh (params from set-params) ##
########################################

# set mesh params, tuned by agg.size
mesh.params <- data.table(size = c(25, 50, 75, 100),
                          max.edge.inner.mult = c(.1, .15, .35, .5))
max.edge = diff(range(boundary.coords[,1]))/(3*5)
bound.outer = diff(range(boundary.coords[,1]))/3

mesh <- fm_mesh_2d(
  #sub.dat[feat == 'Igf2', .(x, y)],
  boundary = list(fm_segm(c(domain))),
  max.edge = c(mesh.params[size == agg.size,
                           max.edge.inner.mult], 2) * max.edge,
  offset = c(max.edge, bound.outer / 1.5),
  cutoff = max.edge / 10,
  min.angle = 21)
#plot(mesh)
# check which mesh vertices are inside the domain
ds1 <- splancs::as.points(boundary.coords)
#  ds1.poly <- ds1[chull(ds1),]
in.dom <- which(splancs::inout(mesh$loc[, 1:2], ds1))

########################################
## run model across selected features ##
########################################

for(lr.n in f.to.mod){

  cat('\n\n\n')
  for(i in 1:3){print(glue('ON FEAT: {lr.n}: {which(f.to.mod==lr.n)} of {length(f.to.mod)}\n'))}
  cat('\n\n\n')

  # make data for this run
  run.dat <- sub.dat[feat == lr.n, ]
  set(run.dat, j = "feat.present", value = as.integer(run.dat[, feat.present]))
  set(run.dat, j = "total.present", value = as.integer(run.dat[, total.present]))

  if(FALSE){
    # visualize raw data
    fields::quilt.plot(run.dat[, x], run.dat[, y], run.dat[, feat.count],
                       nx = run.dat[, length(unique(x))], ny = run.dat[, length(unique(y))])
  }

  # run the model
  source(file.path(c.d, "02-array-run-model.R"))
}
