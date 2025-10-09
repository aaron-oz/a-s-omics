## for a selected feature, this script performs the main spatial
## modeling, generates predictions, calculates performance metrics and
## makes some plots
##
## in particular, this runs a joint spatial model, jointly estimating
## the total feature counts and uses those estimated totals to
## normalize the counts of a particular feature


##########################################
## set params that can't be set earlier ##
##########################################

# predict resolution (number of pixels in x and y)
qp.res.x <- run.dat[, length(unique(x))]
qp.res.y <- run.dat[, length(unique(y))]

dir.create(file.path(o.d, "prediction-objects")) # dir for model fit objects
## dir.create(file.path(o.d, "fitted-models")) # dir for prediction outputs


############################
## setup spde matern and bru model components
############################

## prior matern params
## c(a, b, c, d), where
## P(sp.range < a) = b
## P(sp.sigma > c) = d
matern.pri.total <- c(500, .5, 3, .01) ## a, b, c, d
matern.total <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                                    constr = TRUE, # integrate-to-zero constraint
                                    prior.range = matern.pri.total[1:2],
                                    prior.sigma = matern.pri.total[3:4])

matern.pri.feat.count <- c(500, .5, 3, .01) ## a, b, c, d
matern.feat <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                                   constr = TRUE, # integrate-to-zero constraint
                                   prior.range = matern.pri.feat.count[1:2],
                                   prior.sigma = matern.pri.feat.count[3:4])

comps <- ~ total.int(1) + feat.int(1) +
  total.field(cbind(x, y), model = matern.total) +
  feat.field(cbind(x, y), model = matern.feat)

#########################
## run INLA bru model(s)
#########################

n.samp <- 1000

##~~~~~~~~~~~~~~~~~
## poisson model for feature only
##~~~~~~~~~~~~~~~
## poi.f.lik <- bru_obs(
##   family = "poisson",
##   data = run.dat,
##   formula = feat.count ~ feat.int + feat.field,
##   E = total.count
## )

## fit.f.poi <- bru(
##   comps,
##   poi.f.lik,
##   options = list(bru_verbose = 4)
## )
## summary(fit.f.poi)

## pred.f.poi <- predict(
##   fit.f.poi, run.dat,
##   ~ {
##     feat.density.per.count <- exp( feat.int + feat.field )
##     feat <- total.count * feat.density.per.count
##     list(
##       total = total.count, # data
##       feat.density.per.count = feat.density.per.count,
##       feat = feat
##     )
##   },
##   n.samples = n.samp
## )


##~~~~~~~~~~~~~~~~~
## poisson model for total only
##~~~~~~~~~~~~~~~
## poi.t.lik <- bru_obs(
##   family = "poisson",
##   data = run.dat,
##   formula = total.count ~ total.int + total.field,
##   E = 1
## )

## fit.t.poi <- bru(
##   comps,
##   poi.t.lik,
##   options = list(bru_verbose = 4)
## )
## summary(fit.t.poi)

## pred.t.poi <- predict(
##   fit.t.poi, run.dat,
##   ~ {
##     total <- exp( total.int + total.field )
##     feat.density.per.count <- feat.count / total.count # data
##     feat <- feat.count # data
##     list(
##       total = total,
##       feat.density.per.count = feat.density.per.count,
##       feat = feat
##     )
##   },
##   n.samples = n.samp
## )


##~~~~~~~~~~~~~~~~~
## poisson-poisson joint model
##~~~~~~~~~~~~~~~~~
poi.jt.lik <- inlabru::bru_obs(
  family = "poisson",
  data = run.dat,
  formula = total.count ~ total.int + total.field,
  E = 1
)

poi.jf.lik <- bru_obs(
  family = "poisson",
  data = run.dat,
  formula = feat.count ~ {
    eta.feat.rate.per.count <- feat.int + feat.field
    eta.E.count   <- total.int + total.field
    ## note, formula is expected on link scale - log for poisson, so:
    #      feat.rate = feat.rate.per.count * eta.E.count, and thus
    #      eta.feat = log(feat.rate) = log(feat.rate.per.count) + log(count)
    eta.feat <- eta.feat.rate.per.count + eta.E.count
    eta.feat
  },
  E = 1
)

## run inlabru to approx the posterior
fit.poi.poi <- bru(
  comps,
  poi.jt.lik, poi.jf.lik,
  options = list(
    bru_verbose = 4,
    bru_max_iter = 1,
    control.inla = list(
      int.strategy = "eb"
    )
  )
)
summary(fit.poi.poi)

# sample from the posterior
pred.poi.poi <- predict(
  fit.poi.poi, run.dat,
  ~ {
    total <- exp( total.int + total.field )
    feat.density.per.count <- exp( feat.int + feat.field )
    feat <- total * feat.density.per.count
    list(
      total = total,
      feat.density.per.count = feat.density.per.count,
      feat = feat
    )
  },
  n.samples = n.samp
)

#######################
## make plots
#######################

png(file.path(o.d, glue('data-vs-model-poi-poi-{lr.n}.png')),
    width = (qp.res.x / qp.res.y) * 13 + 8, height = 13, units = 'in', res = 300)
par(mfrow = c(3, 3),
    mai = c(.62, 0.82, .62, 1.22))
fields.style();quilt.plot(run.dat[, x], run.dat[, y], run.dat[, feat.count],
                          main = glue('Observed {lr.n} Counts'),
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x, cex = 1.4)
fields.style();quilt.plot(run.dat[, x], run.dat[, y], run.dat[, feat.count / total.count],
                          main = glue('Observed {lr.n} Counts per Total Counts'),
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
fields.style();quilt.plot(run.dat[, x], run.dat[, y], run.dat[, total.count],
                          main = 'Total Counts',
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
fields.style();quilt.plot(pred.poi.poi$feat[, x],
                          pred.poi.poi$feat[, y],
                          pred.poi.poi$feat[, mean],
                          main = glue('Estimated {lr.n} Counts' ),
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
fields.style();quilt.plot(pred.poi.poi$feat[, x],
                          pred.poi.poi$feat[, y],
                          pred.poi.poi$feat.density.per.count[, mean],
                          main = glue('Estimated {lr.n} Counts per Total Counts' ),
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
fields.style();quilt.plot(pred.poi.poi$total[, x],
                          pred.poi.poi$total[, y],
                          pred.poi.poi$total[, mean],
                          main = glue('Estimated Total Counts'),
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
if(!('feat.count' %in% colnames(pred.poi.poi$feat))){
  pred.poi.poi$feat <-
    merge(pred.poi.poi$feat,
          run.dat[, .(x, y, feat.count, total.count)])
}
fields.style();quilt.plot(pred.poi.poi$feat[, x],
                          pred.poi.poi$feat[, y],
                          pred.poi.poi$feat[, feat.count] - pred.poi.poi$feat[, mean],
                          main = glue('Residual {lr.n} Counts' ),
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
fields.style();quilt.plot(pred.poi.poi$feat[, x],
                          pred.poi.poi$feat[, y],
                          pred.poi.poi$feat[, feat.count / total.count] -
                            pred.poi.poi$feat.density.per.count[, mean],
                          main = glue('Residual {lr.n} Counts per Total Counts' ),
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
fields.style();quilt.plot(pred.poi.poi$total[, x],
                          pred.poi.poi$total[, y],
                          pred.poi.poi$feat[, total.count] - pred.poi.poi$total[, mean],
                          main = glue('Residual Total Counts' ),
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
dev.off()


##############################
## save and cleanup
##############################

cat('\n')
for(i in 1){print(glue('{lr.n}:  saving outputs'))}
cat('\n')

## saving each fit objects takes up about 2gb!
## saveRDS(fit.pois.pois, file = file.path(o.d, "fitted-models",
##                                         glue('fit-poi-poi-{lr.n}.rds')))

saveRDS(pred.poi.poi, file = file.path(o.d, "prediction-objects",
                                       glue('{lr.n}-pred-poi.rds')))


# clean up
rm(run.dat, fit.poi.poi, pred.poi.poi);for(i in 1:3){gc()}
