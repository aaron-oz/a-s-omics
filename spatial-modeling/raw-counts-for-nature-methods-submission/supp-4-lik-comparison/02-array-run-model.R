## for a selected feature, this script performs the main spatial
## modeling, generates predictions, calculates performance metrics and
## makes some plots
##
## in particular, this runs spatial models for estimating feature counts using:
## 1) a poisson likelihood
## 2) a zero-inflated poisson likelihood
## 3) a zero-adjusted poisson likelihood
## and compares the results from all  3

##########################################
## set params that can't be set earlier ##
##########################################

# predict resolution (number of pixels in x and y)
qp.res.x <- run.dat[, length(unique(x))]
qp.res.y <- run.dat[, length(unique(y))]

dir.create(file.path(o.d, "prediction-objects")) # dir for model fit objects
## dir.create(file.path(o.d, "fitted-models")) # dir for prediction outputs


############################
## setup spde matern
############################

## prior matern params
## c(a, b, c, d), where
## P(sp.range < a) = b
## P(sp.sigma > c) = d
#matern.pri.total <- c(2000, .95, 1, .05) ## a, b, c, d
matern.pri.feat.count <- c(500, .95, 1, .05) ## a, b, c, d
matern.pri.feat.present <- c(1000, .95, 10, .05) ## a, b, c, d

## matern.remain <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
##                                      constr = TRUE, # integrate-to-zero constraint
##                                      prior.range = matern.pri.total[1:2],
##                                      prior.sigma = matern.pri.total[3:4])
matern.feat.count <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                                         constr = TRUE, # integrate-to-zero constraint
                                         prior.range = matern.pri.feat.count[1:2],
                                         prior.sigma = matern.pri.feat.count[3:4])
matern.feat.present <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                                           constr = TRUE, # integrate-to-zero constraint
                                           prior.range = matern.pri.feat.present[1:2],
                                           prior.sigma = matern.pri.feat.present[3:4])


#########################
## run INLA bru model(s)
#########################

# We want to obtain CPO data from the estimations
bru_options_set(control.compute = list(cpo = TRUE))
n.samp <- 1000


##~~~~~~~~~~~~~~~~~
## poisson model
##~~~~~~~~~~~~~~~~~
if(poi.mod){

  matern.pri.feat.count <- c(500, .95, .1, .05) ## a, b, c, d
  matern.feat.count <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                                           constr = TRUE, # integrate-to-zero constraint
                                           prior.range = matern.pri.feat.count[1:2],
                                           prior.sigma = matern.pri.feat.count[3:4])

  comps <- ~
    feat.count.int(1) +
    feat.count.field(cbind(x, y), model = matern.feat.count)

  poi.lik <- bru_obs(
    family = "poisson",
    data = run.dat,
    formula = feat.count ~ feat.count.int + feat.count.field,
    E = total.count
  )

  try(
  {
    cat('\n')
    for(i in 1){print(glue('{lr.n}: fitting poisson model'))}
    cat('\n')

    fit.poi <- bru(comps,
                   poi.lik,
                   # options
                   options = list(bru_verbose = 4,
                                  bru_max_iter = 1)
                   )

  }, silent = T)

  if(exists("fit.poi")){
    try(
    {
      cat('\n')
      for(i in 1){print(glue('{lr.n}: predicting from fitted poisson model'))}
      cat('\n')

      pred.poi <- predict(
        fit.poi, run.dat,
        ~ {
          lambda <- exp( feat.count.int + feat.count.field )
          expect <- lambda * total.count
          list(
            lambda = lambda,
            expect = expect,
            obs_prob = dpois(feat.count, expect)
          )
        },
        n.samples = n.samp
      )
    }, silent = T)
  }

  if(exists("fit.poi")){
    poi_pit <- fit.poi$cpo$pit * c(NA_real_, 1)[1 + (run.dat$feat.count > 0)]
  }else{
    poi_pit <- NA
  }

  if(exists("pred.poi")){
    # For Poisson, the posterior conditional variance is equal to
    # the posterior conditional mean, so no need to compute it separately.
    expect_poi <- pred.poi$expect
    expect_poi$pred_var <- expect_poi$mean + expect_poi$sd^2
    expect_poi$log_score <- -log(pred.poi$obs_prob$mean)
  }else{
    expect_poi$expect <- run.dat
    expect_poi$mean <- NA
    expect_poi$median <- NA
    expect_poi$pred_var <- NA
    expect_poi$log_score <- NA
  }
}

##~~~~~~~~~~~~~~~~~
## zip model
##~~~~~~~~~~~~~~~~~

if(zip.mod){

  matern.pri.feat.count <- c(500, .95, .1, .05) ## a, b, c, d
  matern.feat.count <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                                           constr = TRUE, # integrate-to-zero constraint
                                           prior.range = matern.pri.feat.count[1:2],
                                           prior.sigma = matern.pri.feat.count[3:4])

  comps <- ~
    feat.count.int(1) +
    feat.count.field(cbind(x, y), model = matern.feat.count)

  zip.lik <- bru_obs(
    family = "zeroinflatedpoisson1",
    data = run.dat,
    formula = feat.count ~ feat.count.int + feat.count.field,
    E = total.count
  )

  fit.zip <- try(
  {
    cat('\n')
    for(i in 1){print(glue('{lr.n}: fitting zip model'))}
    cat('\n')

    bru(comps,
        zip.lik,
        # options
        options = list(bru_verbose = 4,
                       bru_max_iter = 1)
        )

  }, silent = T)

  if(exists("fit.zip")){
    try(
    {
      cat('\n')
      for(i in 1){print(glue('{lr.n}:  predicting from fitted zip model'))}
      cat('\n')

      pred.zip <- predict(
        fit.zip, run.dat,
        ~ {
          scaling_prob <- (1 - zero_probability_parameter_for_zero_inflated_poisson_1)
          lambda <- exp( feat.count.int + feat.count.field )
          expect_param <- lambda * total.count
          expect <- scaling_prob * expect_param
          variance <- scaling_prob * expect_param *
            (1 + (1 - scaling_prob) * expect_param)
          list(
            lambda = lambda,
            expect = expect,
            variance = variance,
            obs_prob = (1 - scaling_prob) * (feat.count == 0) +
              scaling_prob * dpois(feat.count, expect_param)
          )
        },
        n.samples = n.samp
      )
    }, silent = T)
  }

  if(exists("fit.zip")){
    zip_pit <- fit.zip$cpo$pit * c(NA_real_, 1)[1 + (run.dat$feat.count > 0)]
  }else{
    zip_pit <- NA
  }

  if(exists("pred.zip")){
    expect_zip <- pred.zip$expect
    expect_zip$pred_var <- pred.zip$variance$mean + expect_zip$sd^2
    expect_zip$log_score <- -log(pred.zip$obs_prob$mean)
  }else{
    expect_zip$expect <- run.dat
    expect_zip$mean <- NA
    expect_zip$median <- NA
    expect_zip$pred_var <- NA
    expect_zip$log_score <- NA
  }
}

##~~~~~~~~~~~~~~~~~
## zap model
##~~~~~~~~~~~~~~~~~
if(zap.mod){

  matern.pri.feat.count <- c(500, .95, .1, .05) ## a, b, c, d
  matern.pri.feat.present <- c(1000, .95, 10, .05) ## a, b, c, d
  matern.feat.count <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                                           constr = TRUE, # integrate-to-zero constraint
                                           prior.range = matern.pri.feat.count[1:2],
                                           prior.sigma = matern.pri.feat.count[3:4])
  matern.feat.present <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                                             constr = TRUE, # integrate-to-zero constraint
                                             prior.range = matern.pri.feat.present[1:2],
                                             prior.sigma = matern.pri.feat.present[3:4])


  comps <- ~
    ## for modeling the total
    ##  total.count.int(1) +
    ##  total.count.field(cbind(x, y), model = matern.feat.present) +
    # for modeling feature presence/absence
    feat.present.int(1) +
    feat.present.field(cbind(x, y), model = matern.feat.present) +
    # for modeling feature counts (when present is using a zip/zap)
    feat.count.int(1) +
    feat.count.field(cbind(x, y), model = matern.feat.count)

  zap.feat.lik <- bru_obs(
    family = "nzpoisson",
    data = run.dat[feat.count > 0, ],
    # feat.count.lambda  = exp(feat.count.eta)*total.count.lambda
    # feat.count.lambda  = exp(feat.count.eta)*exp(total.count.eta)
    # feat.count.lambda = exp(feat.count.eta + total.count.eta)
    formula = feat.count ~ feat.count.int + feat.count.field, #+ total.count.int + total.count.field,
    E = total.count
  )

  zap.pres.lik <- bru_obs(
    family = "binomial",
    data = run.dat,
    formula = feat.present ~ feat.present.int + feat.present.field
  )

  try(
  {
    cat('\n')
    for(i in 1){print(glue('{lr.n}: fitting zap model'))}
    cat('\n')

    fit.zap <- bru(comps,
                   zap.feat.lik,
                   zap.pres.lik,
                   # options
                   options = list(bru_verbose = 4,
                                  bru_max_iter = 1)
                   )
  }, silent = T)

  if(exists("fit.zap")){
    try(
    {
      cat('\n')
      for(i in 1){print(glue('{lr.n}: predicting from fitted zap model'))}
      cat('\n')

      pred.zap <- predict(
        fit.zap,
        run.dat,
        ~ {
          presence_prob <- plogis( feat.present.int + feat.present.field )
          lambda <- exp(feat.count.int + feat.count.field)
          expect_param <- presence_prob * lambda * total.count
          expect <- expect_param / (1 - exp(-lambda * total.count))
          variance <- expect * (1 - exp(-lambda * total.count) * expect)
          list(
            presence = presence_prob,
            lambda = lambda,
            expect = expect,
            variance = variance,
            obs_prob = (1 - presence_prob) * (feat.count == 0) +
              (feat.count > 0) * presence_prob * dpois(feat.count, expect_param) /
              (1 - dpois(0, expect_param))
          )
        },
        n.samples = n.samp
      )
    }, silent = T)
  }

  if(exists("fit.zap")){
    zap_pit <- rep(NA_real_, nrow(run.dat))
    zap_pit[run.dat$feat.count > 0] <- fit.zap$cpo$pit[-seq_len(nrow(run.dat))]
  }

  if(exists("pred.zap")){
    presence_zap <- pred.zap$presence
    expect_zap <- pred.zap$expect
    expect_zap$pred_var <- pred.zap$variance$mean + expect_zap$sd^2
    expect_zap$log_score <- -log(pred.zap$obs_prob$mean)
  }else{
    expect_zap$expect <- run.dat
    expect_zap$mean <- NA
    expect_zap$median <- NA
    expect_zap$pred_var <- NA
    expect_zap$log_score <- NA
  }
}

######################
## model comparison
######################
if(run.lik.comp){
  cat('\n')
  for(i in 1){print(glue('{lr.n}:  model comparison'))}
  cat('\n')

  comp.df <- data.frame(
    count = rep(run.dat$feat.count, times = 3),
    pred_mean = c(
      expect_poi$mean,
      expect_zip$mean,
      expect_zap$mean
    ),
    pred_var = c(
      expect_poi$pred_var,
      expect_zip$pred_var,
      expect_zap$pred_var
    ),
    pred_median = c(
      expect_poi$median,
      expect_zip$median,
      expect_zap$median
    ),
    log_score = c(
      expect_poi$log_score,
      expect_zip$log_score,
      expect_zap$log_score
    ),
    pit = c(
      poi_pit,
      zip_pit,
      zap_pit
    ),
    Model = rep(c("Poisson", "ZIP", "ZAP"), each = nrow(run.dat))
  )

  resid.p <- ggplot(comp.df) +
    geom_point(aes(pred_mean, count - pred_mean, color = Model)) +
    ggtitle("Residuals")

  pit.p <- ggplot(comp.df) +
    stat_ecdf(aes(pit, color = Model), na.rm = TRUE) +
    scale_x_continuous(expand = c(0, 0)) +
    ggtitle("PIT")

  diag.p <- patchwork::wrap_plots(resid.p, pit.p, nrow = 1, guides = "collect")

  comp.df <- comp.df %>%
    mutate(
      AE = abs(count - pred_median),
      SE = (count - pred_mean)^2,
      DS = (count - pred_mean)^2 / pred_var + log(pred_var),
      LG = log_score
    )

  scores <- comp.df %>%
    group_by(Model) %>%
    summarise(
      MAE = mean(AE, na.rm = T),
      RMSE = sqrt(mean(SE, na.rm = T)),
      MDS = mean(DS, na.rm = T),
      MLG = mean(LG, na.rm = T)
    ) %>%
    left_join(
      data.frame(
        Model = c("Poisson", "ZIP", "ZAP"),
        Order = 1:3
      ),
      by = "Model"
    ) %>%
    arrange(Order) %>%
    select(-Order)
  knitr::kable(scores)


  #######################
  ## make plots
  #######################

  ## TODO start here, make comparison plots

  # poi (-, density, count, count resid)
  # zip (-, density, count, count resid)
  # zap (present, density, count, count resid)
  # obs (present, density, count)

  ## TODO then, update save outputs section, and run

  cat('\n')
  for(i in 1){print(glue('{lr.n}:  making plots'))}
  cat('\n')

  ## first make plots with constant color range per column

  prob.zlim <- c(0, 1)
  prob.cols <- (magma(256))

  dens.zlim <- c(ifelse(exists("pred.poi"), NA, pred.poi$lambda$mean),
                 ifelse(exists("pred.zip"), NA, pred.zip$lambda$mean),
                 ifelse(exists("pred.zap"), NA, pred.zap$lambda$mean),
                 run.dat[, feat.count / total.count]) |> range(na.rm = T)
  dens.cols <- cividis(256)

  cnt.zlim <- c(ifelse(exists("pred.poi"), NA, pred.poi$expect$mean),
                ifelse(exists("pred.zip"), NA, pred.zip$expect$mean),
                ifelse(exists("pred.zap"), NA, pred.zap$expect$mean),
                run.dat[, feat.count]) |> range(na.rm = T)
  cnt.cols <- viridis(256)

  res.zlim <- c(ifelse(exists("pred.poi"), NA, run.dat$feat.count - pred.poi$expect$mean),
                ifelse(exists("pred.zip"), NA, run.dat$feat.count - pred.zip$expect$mean),
                ifelse(exists("pred.zap"), NA, run.dat$feat.count - pred.zap$expect$mean),
                run.dat[, feat.count / total.count]) |> range(na.rm = T)
  res.cols <- turbo(256)

  png(file.path(o.d, glue('{lr.n}-model-prediction-comparisons-fixed-colors.png')),
      width = (qp.res.x / qp.res.y) * 4 / 3 * 13 + 8, height = 4 / 3 * 13, units = 'in', res = 300)
  par(mfrow = c(4, 4),
      mai = c(.62, 0.82, .62, 1.22))

  if(exists("pred.poi")){

    fields.style();quilt.plot(pred.poi$expect$x,
                              pred.poi$expect$y,
                              1 - dpois(0, pred.poi$expect$mean),
                              main = glue('{lr.n} - Poisson - Estimated Prob of Occurence'),
                              zlim = prob.zlim, nlevel = 256, col = prob.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.poi$lambda$x,
                              pred.poi$lambda$y,
                              pred.poi$lambda$mean,
                              main = glue('{lr.n} - Poisson - Estimated Density'),
                              zlim = dens.zlim, nlevel = 256, col = dens.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.poi$expect$x,
                              pred.poi$expect$y,
                              pred.poi$expect$mean,
                              main = glue('{lr.n} - Poisson - Estimated Count'),
                              zlim = cnt.zlim, nlevel = 256, col = cnt.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.poi$expect$x,
                              pred.poi$expect$y,
                              run.dat$feat.count - pred.poi$expect$mean,
                              main = glue('{lr.n} - Poisson - Count Residuals'),
                              zlim = res.zlim, nlevel = 256, col = res.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  }else{
    for(i in 1:4){ plot.new() }
  }

  if(exists("pred.zip")){
    fields.style();quilt.plot(pred.zip$expect$x,
                              pred.zip$expect$y,
                              1 - dpois(0, pred.zip$expect$mean),
                              main = glue('{lr.n} - ZIP - Estimated Prob of Occurence'),
                              zlim = prob.zlim, nlevel = 256, col = prob.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.zip$lambda$x,
                              pred.zip$lambda$y,
                              pred.zip$lambda$mean,
                              main = glue('{lr.n} - ZIP - Estimated Density'),
                              zlim = dens.zlim, nlevel = 256, col = dens.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.zip$expect$x,
                              pred.zip$expect$y,
                              pred.zip$expect$mean,
                              main = glue('{lr.n} - ZIP - Estimated Count'),
                              zlim = cnt.zlim, nlevel = 256, col = cnt.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.zip$expect$x,
                              pred.zip$expect$y,
                              run.dat$feat.count - pred.zip$expect$mean,
                              main = glue('{lr.n} - ZIP - Count Residuals'),
                              zlim = res.zlim, nlevel = 256, col = res.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  }else{
    for(i in 1:4){ plot.new() }
  }

  if(exists("pred.zap")){
    fields.style();quilt.plot(pred.zap$expect$x,
                              pred.zap$expect$y,
                              1 - dpois(0, pred.zap$expect$mean),
                              main = glue('{lr.n} - ZAP - Estimated Prob of Occurence'),
                              zlim = prob.zlim, nlevel = 256, col = prob.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.zap$lambda$x,
                              pred.zap$lambda$y,
                              pred.zap$lambda$mean,
                              main = glue('{lr.n} - ZAP - Estimated Density'),
                              zlim = dens.zlim, nlevel = 256, col = dens.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.zap$expect$x,
                              pred.zap$expect$y,
                              pred.zap$expect$mean,
                              main = glue('{lr.n} - ZAP - Estimated Count'),
                              zlim = cnt.zlim, nlevel = 256, col = cnt.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.zap$expect$x,
                              pred.zap$expect$y,
                              run.dat$feat.count - pred.zap$expect$mean,
                              main = glue('{lr.n} - ZAP - Count Residuals'),
                              zlim = res.zlim, nlevel = 256, col = res.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  }else{
    for(i in 1:4){ plot.new() }
  }

  fields.style();quilt.plot(run.dat$x,
                            run.dat$y,
                            run.dat$feat.present,
                            main = glue('{lr.n} - Data - Prob of Occurence'),
                            zlim = prob.zlim, nlevel = 256, col = prob.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(run.dat$x,
                            run.dat$y,
                            run.dat[, feat.count / total.count],
                            main = glue('{lr.n} - Data -  Density'),
                            zlim = dens.zlim, nlevel = 256, col = dens.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(run.dat$x,
                            run.dat$y,
                            run.dat[, feat.count],
                            main = glue('{lr.n} - Data -  Counts'),
                            zlim = cnt.zlim, nlevel = 256, col = cnt.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)

  fields.style();quilt.plot(run.dat$x,
                            run.dat$y,
                            run.dat[, total.count],
                            main = 'Data -  NUMI',
                            nlevel = 256, col = res.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  dev.off()


  png(file.path(o.d, glue('{lr.n}-model-prediction-comparisons-free-colors.png')),
      width = (qp.res.x / qp.res.y) * 4 / 3 * 13 + 8, height = 4 / 3 * 13, units = 'in', res = 300)
  par(mfrow = c(4, 4),
      mai = c(.62, 0.82, .62, 1.22))

  if(exists("pred.poi")){
    fields.style();quilt.plot(pred.poi$expect$x,
                              pred.poi$expect$y,
                              1 - dpois(0, pred.poi$expect$mean),
                              main = glue('{lr.n} - Poisson - Estimated Prob of Occurence'),
                              nlevel = 256, col = prob.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.poi$lambda$x,
                              pred.poi$lambda$y,
                              pred.poi$lambda$mean,
                              main = glue('{lr.n} - Poisson - Estimated Density'),
                              nlevel = 256, col = dens.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.poi$expect$x,
                              pred.poi$expect$y,
                              pred.poi$expect$mean,
                              main = glue('{lr.n} - Poisson - Estimated Count'),
                              nlevel = 256, col = cnt.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.poi$expect$x,
                              pred.poi$expect$y,
                              run.dat$feat.count - pred.poi$expect$mean,
                              main = glue('{lr.n} - Poisson - Count Residuals'),
                              nlevel = 256, col = res.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  }else{
    for(i in 1:4){ plot.new() }
  }

  if(exists("pred.zip")){
    fields.style();quilt.plot(pred.zip$expect$x,
                              pred.zip$expect$y,
                              1 - dpois(0, pred.zip$expect$mean),
                              main = glue('{lr.n} - ZIP - Estimated Prob of Occurence'),
                              nlevel = 256, col = prob.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.zip$lambda$x,
                              pred.zip$lambda$y,
                              pred.zip$lambda$mean,
                              main = glue('{lr.n} - ZIP - Estimated Density'),
                              nlevel = 256, col = dens.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.zip$expect$x,
                              pred.zip$expect$y,
                              pred.zip$expect$mean,
                              main = glue('{lr.n} - ZIP - Estimated Count'),
                              nlevel = 256, col = cnt.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.zip$expect$x,
                              pred.zip$expect$y,
                              run.dat$feat.count - pred.zip$expect$mean,
                              main = glue('{lr.n} - ZIP - Count Residuals'),
                              nlevel = 256, col = res.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  }else{
    for(i in 1:4){ plot.new() }
  }

  if(exists("pred.zap")){
    fields.style();quilt.plot(pred.zap$expect$x,
                              pred.zap$expect$y,
                              1 - dpois(0, pred.zap$expect$mean),
                              main = glue('{lr.n} - ZAP - Estimated Prob of Occurence'),
                              nlevel = 256, col = prob.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.zap$lambda$x,
                              pred.zap$lambda$y,
                              pred.zap$lambda$mean,
                              main = glue('{lr.n} - ZAP - Estimated Density'),
                              nlevel = 256, col = dens.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.zap$expect$x,
                              pred.zap$expect$y,
                              pred.zap$expect$mean,
                              main = glue('{lr.n} - ZAP - Estimated Count'),
                              nlevel = 256, col = cnt.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.zap$expect$x,
                              pred.zap$expect$y,
                              run.dat$feat.count - pred.zap$expect$mean,
                              main = glue('{lr.n} - ZAP - Count Residuals'),
                              nlevel = 256, col = res.cols,
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  }else{
    for(i in 1:4){ plot.new() }
  }

  fields.style();quilt.plot(run.dat$x,
                            run.dat$y,
                            run.dat$feat.present,
                            main = glue('{lr.n} - Data - Prob of Occurence'),
                            nlevel = 256, col = prob.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(run.dat$x,
                            run.dat$y,
                            run.dat[, feat.count / total.count],
                            main = glue('{lr.n} - Data -  Density'),
                            nlevel = 256, col = dens.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(run.dat$x,
                            run.dat$y,
                            run.dat[, feat.count],
                            main = glue('{lr.n} - Data -  Counts'),
                            nlevel = 256, col = cnt.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(run.dat$x,
                            run.dat$y,
                            run.dat[, total.count],
                            main = 'Data -  NUMI',
                            nlevel = 256, col = res.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  dev.off()



  ## poi.p.d <-  ggplot() +
  ##   geom_fm(data = mesh) +
  ##   geom_raster(data = expect_poi, aes(x, y, fill = mean / total.count)) +
  ##   ggtitle("Poisson predicted density")

  ## poi.p.c <-  ggplot() +
  ##   geom_fm(data = mesh) +
  ##   geom_raster(data = expect_poi, aes(x, y, fill = mean)) +
  ##   ggtitle("Poisson predicted counts")

  ## zip.p.d <-  ggplot() +
  ##   geom_fm(data = mesh) +
  ##   geom_raster(data = expect_zip, aes(x, y, fill = mean / total.count)) +
  ##   ggtitle("ZIP predicted density")

  ## zip.p.c <-  ggplot() +
  ##   geom_fm(data = mesh) +
  ##   geom_raster(data = expect_zip, aes(x, y, fill = mean)) +
  ##   ggtitle("ZIP predicted counts")

  ## zap.p.p <- ggplot() +
  ##   geom_fm(data = mesh) +
  ##   geom_raster(data = presence_zap, aes(x, y, fill = mean)) +
  ##   geom_point(data = run.dat[feat.count > 0, ], aes(x, y),
  ##              color = "firebrick", size = 1, pch = 4, alpha = 1) +
  ##   ggtitle("Presence probability")
  ## zap.p.d <- ggplot() +
  ##   geom_fm(data = mesh) +
  ##   geom_raster(data = expect_zap, aes(x, y, fill = mean / total.count)) +
  ##   ggtitle("Expected feature counts per total count")
  ## # geom_sf(data = nests, color = "firebrick", size = 1, pch = 4, alpha = 0.2) +
  ## zap.p.c <- ggplot() +
  ##   geom_fm(data = mesh) +
  ##   geom_raster(data = expect_zap, aes(x, y, fill = mean)) +
  ##   #    geom_sf(data = nests, color = "firebrick", size = 1, pch = 4, alpha = 0.2) +
  ##   ggtitle("Expected feature counts per total count")

  ## patchwork::wrap_plots(zap.p.p, zap.p.d, zap.p.c, nrow = 1)
}





cat('\n')
for(i in 1){print(glue('{lr.n}:  saving outputs'))}
cat('\n')

## saving each fit objects takes up about 2gb!
## saveRDS(fit.pois.pois, file = file.path(o.d, "fitted-models",
##                                         glue('fit-pois-pois-{lr.n}.rds')))

# save model diagnostic plots
png(file.path(o.d, glue("{lr.n}-diagnostic-plots.png")), width = 10, height = 4, units = 'in', res = 100)
print(diag.p)
dev.off()

# save model comparison metrics
fwrite(scores, file = file.path(o.d, glue("{lr.n}-scores.csv")))

if(exists("pred.poi")){
  saveRDS(pred.poi, file = file.path(o.d, "prediction-objects",
                                     glue('{lr.n}-pred-poi.rds')))
}

if(exists("pred.zip")){
  saveRDS(pred.zip, file = file.path(o.d, "prediction-objects",
                                     glue('{lr.n}-pred-zip.rds')))
}

if(exists("pred.zap")){
  saveRDS(pred.zap, file = file.path(o.d, "prediction-objects",
                                     glue('{lr.n}-pred-zap.rds')))
}

# clean up
rm(run.dat, fit.poi, fit.zip, fit.zap,
   pred.poi, pred.zip, pred.zap,
   expect_poi, expect_zip, expect_zap);gc()
