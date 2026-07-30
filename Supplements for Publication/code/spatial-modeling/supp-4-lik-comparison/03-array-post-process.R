## this script is meant to be run after all features have been modeled
## and is used to pull together results from across the independent
## feature runs.

## setwd("/path/to/repo/")
## i.d <- file.path(getwd(), "/data-inputs/mouse-embryo-raw")
## o.d <- "/path/to/repo/data-outputs/mouse-embryo-raw/run-date" # initial output location
## require(data.table)
## require(glue)
## require(fields)

# make plots?
plot.preds <- T

# get files in output dir
a.f.p <- list.files(file.path(o.d, "prediction-objects"), full.names = T)
a.f   <- list.files(file.path(o.d, "prediction-objects"))

# subset to poisson pred objects
p.f.p <- grep("poi.rds", a.f.p, value = T)
p.f   <- grep("poi.rds", a.f, value = T)

# get list of all LR mechanism features
lr.feats <- grep(".rds", p.f, value = T) |>
  strsplit("\\.") |>
  lapply(function(x){x[[1]]}) |>
  unlist() |>
  strsplit("-") |>
  lapply(function(x){x[[1]]}) |>
  unlist() |> unique() |> sort()

# and get the individual L and R names
lr.feats <- strsplit(lr.feats, "-") |> unlist() |> unique() |> sort()

# load in LR pairs
load(file.path(i.d, "fantom-hierarchy.Robj"))
hier <- as.data.table(output)
lr.mechs <- data.table(l = hier$LIGAND,
                       r = hier$RECEPTOR,
                       m = hier$MECHANISM)

lr.mechs <- subset(lr.mechs, l %in% toupper(lr.feats) &
                               r %in% toupper(lr.feats))

####################################################
## Collate some objects from across feat results
####################################################

# loop through and get metrics
res.list <- list()
rel.list <- list()
for(ff in lr.feats){
  tmp <- fread(file.path(o.d, glue('{ff}-scores.csv')))
  Feature <- rep(ff, nrow(tmp))
  tmp <- cbind(Feature, tmp)

  rel <- copy(tmp)
  rel[2, 3:6] <- rel[2, 3:6] / rel[1, 3:6]
  rel[3, 3:6] <- rel[3, 3:6] / rel[1, 3:6]

  res.list[[which(lr.feats == ff)]] <- tmp
  rel.list[[which(lr.feats == ff)]] <- rel
}
res.dt <- rbindlist(res.list)
rel.dt <- rbindlist(rel.list)

# summarize results, excluding those with abnormal values indicating
# nonconvergence
sum.res <- res.dt[MAE < 100 & # excludes 4
                    abs(MDS) < 100, # excludes 3
                  .(MAE = mean(MAE),
                    RMSE = mean(RMSE),
                    MDS = mean(MDS),
                    MLG = mean(MLG)),
                  by = Model]

sum.rel <- rel.dt[MAE < 100 & # excludes 4
                    abs(MDS) < 100, # excludes 3
                  .(MAE = mean(MAE),
                    RMSE = mean(RMSE),
                    MDS = mean(MDS),
                    MLG = mean(MLG)),
                  by = Model]

print(xtable::xtable(sum.res))
print(xtable::xtable(sum.rel))

## plot output prediction objects ##
####################################

if(plot.preds){
  # loop through the mechanism preds and process
  dir.create(file.path(o.d, "convolution-plots"))
  dir.create(file.path(o.d, "convolution-data"))

  for(m.idx in 1:lr.mechs[, .N]){

    message(glue("on {m.idx} of {lr.mechs[,.N]}"))

    m.n <- lr.mechs[m.idx, m]
    l.n <- lr.mechs[m.idx, l]
    l.fn <- grep(l.n, p.f.p, ignore.case = T, value = T)
    if(length(l.fn) > 1){
      l.fn <- l.fn[which(unlist(lapply(strsplit(basename(l.fn), "-"), function(x){toupper(x[[1]])})) == l.n)]
    }
    l.pred <- readRDS(l.fn)

    r.n <- lr.mechs[m.idx, r]
    r.fn <- grep(r.n, p.f.p, ignore.case = T, value = T)
    if(length(r.fn) > 1){
      r.fn <- r.fn[which(unlist(lapply(strsplit(basename(r.fn), "-"), function(x){toupper(x[[1]])})) == r.n)]
    }
    r.pred <- readRDS(r.fn)

    calc.conv <- function(l, r,
                          calc.f,
                          k.d = NULL){

      if(calc.f == "min"){
        both <- cbind(as.vector(l), as.vector(r))
        elem.min <- apply(both, 1, min)
        return(elem.min)
      }

      if(calc.f == "geom.mean"){
        return(sqrt(l * r))
      }

      if(calc.f == "prod"){
        return(l * r)
      }

      if(calc.f == "kinetic"){
        if(is.null(k.d) | k.d < 0){
          stop("supply valid, non-neg, k.d for kinetics equilibrium calc")
        }
        # from https://en.wikipedia.org/wiki/Receptor%E2%80%93ligand_kinetics
        E <- (r - l - k.d) / 2
        D <- sqrt(E^2 + r * k.d)
        # take the one positive root to find equilibrium r
        roots <- cbind(E + D, E - D)
        r.e <- apply(roots, 1, max)
        # binding concentration
        c <- r - r.e
        return(c)
      }
    }

    # make one overall results object to use in plotting
    res <- data.table(
      # first, data
      l.n = l.n,
      r.n = r.n,
      m.n = m.n,
      x = l.pred[[1]]$x,
      y = l.pred[[1]]$y,
      t = l.pred[[1]]$total.count,
      l = l.pred[[1]]$feat.count,
      r = r.pred[[1]]$feat.count,
      l.d = l.pred[[1]]$feat.count / l.pred[[1]]$total.count,
      r.d = r.pred[[1]]$feat.count / r.pred[[1]]$total.count,
      # then estimates
      l.est = l.pred$expect$median,
      r.est = r.pred$expect$median,
      # and density estimates
      l.est.d = l.pred$expect$median / l.pred$expect$total.count,
      r.est.d = r.pred$expect$median / r.pred$expect$total.count
    )
    ## and we calculate the convolution many ways
    # first, data
    res[, m.data.min  := calc.conv(l, r, "min")]
    res[, m.data.prod := calc.conv(l, r, "prod")]
    res[, m.data.gm   := calc.conv(l, r, "geom.mean")]
    res[, m.data.k.01 := calc.conv(l, r, "kinetic", k.d = .01)]
    res[, m.data.k.5  := calc.conv(l, r, "kinetic", k.d = .5)]
    res[, m.data.k1   := calc.conv(l, r, "kinetic", k.d = 1)]
    # using data density
    res[, m.data.d.min  := calc.conv(l.d, r.d, "min")]
    res[, m.data.d.prod := calc.conv(l.d, r.d, "prod")]
    res[, m.data.d.gm   := calc.conv(l.d, r.d, "geom.mean")]
    res[, m.data.d.k.01 := calc.conv(l.d, r.d, "kinetic", k.d = .01)]
    res[, m.data.d.k.5  := calc.conv(l.d, r.d, "kinetic", k.d = .5)]
    res[, m.data.d.k1   := calc.conv(l.d, r.d, "kinetic", k.d = 1)]
    # using estimates
    res[, m.est.min  := calc.conv(l.est, r.est, "min")]
    res[, m.est.prod := calc.conv(l.est, r.est, "prod")]
    res[, m.est.gm   := calc.conv(l.est, r.est, "geom.mean")]
    res[, m.est.k.01 := calc.conv(l.est, r.est, "kinetic", k.d = .01)]
    res[, m.est.k.5  := calc.conv(l.est, r.est, "kinetic", k.d = .5)]
    res[, m.est.k1   := calc.conv(l.est, r.est, "kinetic", k.d = 1)]
    # using estimated density
    res[, m.est.d.min  := calc.conv(l.est.d, r.est.d, "min")]
    res[, m.est.d.prod := calc.conv(l.est.d, r.est.d, "prod")]
    res[, m.est.d.gm   := calc.conv(l.est.d, r.est.d, "geom.mean")]
    res[, m.est.d.k.01 := calc.conv(l.est.d, r.est.d, "kinetic", k.d = .01)]
    res[, m.est.d.k.5  := calc.conv(l.est.d, r.est.d, "kinetic", k.d = .5)]
    res[, m.est.d.k1   := calc.conv(l.est.d, r.est.d, "kinetic", k.d = 1)]

    fwrite(res, file = file.path(o.d, "convolution-data", glue("convolution-calcs-{m.n}.csv")))


    ## make plots

    # L   conv1 conv2 conv3
    # R   conv4 conv5 conv6
    # L/T conv1 conv2 conv3
    # R/T conv1 conv2 conv3

    qp.res.x <- res[, length(unique(x))]
    qp.res.y <- res[, length(unique(y))]

    png(file.path(o.d, "convolution-plots", glue('convolution-comparison-poi-{l.n}-{r.n}.png')),
        width = (qp.res.x / qp.res.y) * 13 + 8, height = 20, units = 'in', res = 300)
    par(mfrow = c(6, 4),
        mai = c(.62, 0.82, .62, 1.22))

    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, l.est],
                              main = glue('Estimated {l.n} Counts'),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)
    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, calc.conv(l.est, r.est, "prod")],
                              main = glue('Product of {l.n} and {r.n}'),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)

    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, calc.conv(l.est, r.est, "geom.mean")],
                              main = glue('Geometric Mean of {l.n} and {r.n}'),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)

    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, calc.conv(l.est, r.est, "min")],
                              main = glue('Minimum of {l.n} and {r.n}'),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, r.est],
                              main = glue('Estimated {r.n} Counts'),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)

    for(kk in c(.1, .5, 1)){
      fields.style();quilt.plot(res[, x],
                                res[, y],
                                res[, calc.conv(l.est, r.est, "kinetic", k.d = kk)],
                                main = glue('LR Complex Equil. Conc. of {l.n} and {r.n} with k.d={kk}'),
                                nx = qp.res.x, ny = qp.res.y, asp = 1)
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, l.est.d],
                              main = glue('Estimated {l.n} Density'),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)
    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, calc.conv(l.est.d, r.est.d, "prod")],
                              main = glue('Product of {l.n} and {r.n} Densities'),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)

    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, calc.conv(l.est.d, r.est.d, "geom.mean")],
                              main = glue('Geometric Mean of {l.n} and {r.n} Densities'),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)

    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, calc.conv(l.est.d, r.est.d, "min")],
                              main = glue('Minimum of {l.n} and {r.n} Densities'),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, r.est.d],
                              main = glue('Estimated {r.n} Densities' ),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)

    for(kk in c(.1, .5, 1)){
      fields.style();quilt.plot(res[, x],
                                res[, y],
                                res[, calc.conv(l.est.d, r.est.d, "kinetic", k.d = kk)],
                                main = glue('LR Complex Equil. Conc. of {l.n} and {r.n} Densities with k.d={kk}'),
                                nx = qp.res.x, ny = qp.res.y, asp = 1)
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, l],
                              main = glue('{l.n} Data'),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)
    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, calc.conv(l, r, "prod")],
                              main = glue('Product of {l.n} and {r.n} Data'),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)

    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, calc.conv(l, r, "geom.mean")],
                              main = glue('Geometric Mean of {l.n} and {r.n} Data'),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)

    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, calc.conv(l, r, "min")],
                              main = glue('Minimum of {l.n} and {r.n} Data'),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    fields.style();quilt.plot(res[, x],
                              res[, y],
                              res[, r],
                              main = glue('{r.n} Data' ),
                              nx = qp.res.x, ny = qp.res.y, asp = 1)

    for(kk in c(.1, .5, 1)){
      fields.style();quilt.plot(res[, x],
                                res[, y],
                                res[, calc.conv(l, r, "kinetic", k.d = kk)],
                                main = glue('LR Complex Equil. Conc. of {l.n} and {r.n} Data with k.d={kk}'),
                                nx = qp.res.x, ny = qp.res.y, asp = 1)
    }

    dev.off()
  }
}
