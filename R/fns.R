# %% ####################################################
#' fit CEF using specified algorithm
#' @param estimator string with name of conditional mean estimator
#' @param fmla formula for conditional mean estimator
#' @param df dataframe
#' @param y name of outcome variable (required for packages that don't use formulas to fit)
#' @return outcome model object
#' @export
mhatter = function(estimator, fmla, df, y, ...){
  require(glmnet); require(ranger); require(xgboost)
  # ols
  if(estimator == 'ols'){
    m = lm(fmla, data = df)
  } else if(estimator == 'logit'){
    m = glm(fmla, family = binomial(), data = df)
  }
  else if (estimator %in% c('lasso', 'ridge', 'rforest', 'xgboost') ){
    # models that demand model matrices
    X = model.matrix(fmla, df) # expand out interactions / factors
    ym = df[[y]] |> as.matrix()
    if(estimator == 'lasso'){
      m = cv.glmnet(x = X, y = ym, alpha = 1, ...)
    } else if(estimator == 'ridge'){
      m = cv.glmnet(x = X, y = ym, alpha = 0, ...)
    # random forest
    } else if(estimator == 'rforest'){
      m = ranger(x = X, y = ym)
    } else if(estimator == 'xgboost'){
      m = xgboost(data = X, label = ym, nrounds = 1000, verbose = 0,  # silent,
        # stop if no improvement for 10 consecutive trees)
        early_stopping_rounds = 10)
    }
    #################################################
    # add more here
    #################################################
  }
  # return model object - we will predict from it on entire sample
  return(m)
}

# %% ####################################################
#' fit pscore using specified algorithm
#' @param estimator string with name of pscore estimator
#' @param fmla formula for conditional mean estimator
#' @param df dataframe
#' @param y name of outcome variable
#' @return pscore model object
#' @export
ehatter = function(estimator, fmla, df, w, ...){
  require(glmnet); require(ranger); require(xgboost)
  # logit
  if(estimator == 'logit'){
    m = glm(fmla, family = binomial(), data = df)
  }
    # models that demand model matrices
  if (estimator %in% c('lasso', 'ridge', 'rforest', 'xgboost')){
      X = model.matrix(fmla, df)
      wm = df[[w]]  |> as.matrix()
    if(estimator == 'lasso'){
      m = cv.glmnet(x = X, y = wm, family = 'binomial',
        type.measure = 'class', alpha = 1, ...)
    } else if(estimator == 'ridge'){
      m = cv.glmnet(x = X, y = wm, family = 'binomial',
        type.measure = 'class', alpha = 0, ...)
    } else if(estimator == 'rforest'){
      m = ranger(x = X, y = wm)
    } else if(estimator == 'xgboost'){
      m = xgboost(data = X, label = wm,
        nrounds = 1000, verbose = 0,  # silent,
        # stop if no improvement for 10 consecutive trees)
        early_stopping_rounds = 10)
    }
    #################################################
    # add more here
    #################################################
  }
  return(m)
}

# %% ####################################################
#' estimate ATE by IPW
#' @param est name of estimator
#' @param w name of treatment in dataframe
#' @param y name of outcome
#' @param df dataframe
#' @param fml formula for propensity score
#' @param ret_est boolean for whether to return ATE estimate, returns pscore vectors when false
#' @param psrange range of propensity scores to use for ATE estimation
#' @return vector of treatment effects or ATE estimate
#' @import glmnet
#' @export
ate_ipw = function(est, w, y, df, fml, ret_est = T, psrange = c(0, 1), ...) {
  ehat = ehatter(est, fmla = fml, df = df, w = w, ...)
  if (est %in% c('lasso', 'ridge')){
    # feed predict.glmnet model matrices instead of df predict at lambda.min
    pscore = predict(ehat, model.matrix(fml, df), s = 'lambda.min', type = 'response')
  } else if (est == 'rforest'){
    pscore = predict(ehat, model.matrix(fml, df), type = 'response')[['predictions']]
  } else if (est == 'xgboost'){
    pscore = predict(ehat, model.matrix(fml, df))
  } else {
    pscore = predict(ehat, df, type = 'response')
  }
  # other predict methods go here
  df2     = data.frame(df[[y]], df[[w]], pscore)
  colnames(df2) = c("y", "w", "ehat")
  if(ret_est == T){ # for aipw
    estsamp = df2[df2$ehat >= psrange[1] & df2$ehat <= psrange[2],]
    ate_ipw = with(estsamp, mean((y * w)/ehat - (y * (1-w))/(1-ehat)) )
    return(ate_ipw)
  } else { # for aipw
    return(df2$ehat)
  }
}

# %% ####################################################
#' function to compute ate using regression adjustment
#' @param est string estimator name, must be %in% c("ols", "logit",  "lasso", "ridge", "rforest", "xgboost")
#' @param name of treatment
#' @param name of outcome
#' @param df data.frame
#' @param fml formula for outcome model (cannot include treatment)
#' @param ret_est boolean for whether to return ATE estimate, returns yhat vectors when false
#' @export
ate_reg = function(est, w, y, df, fml, ret_est = T , ...){
  # treatment index
  treatInd = which(df[[w]] == 1)
  # compute CEF model objects
  est_treat = mhatter(est, fml, df[treatInd,  ]  , y = y, ...)
  est_ctrl  = mhatter(est, fml, df[-(treatInd), ], y = y, ...)
  if (est %in% c('lasso', 'ridge')){
    m1 = predict(est_treat, model.matrix(fml, df), s = 'lambda.min', type = 'response')
    m0 = predict(est_ctrl,  model.matrix(fml, df), s = 'lambda.min', type = 'response')
  } else if (est %in% c('rforest', 'xgboost')) {
    m1 = predict(est_treat, model.matrix(fml, df))
    m0 = predict(est_ctrl,  model.matrix(fml, df))
  } else {
    # predict for entire sample
    m1 = predict(est_treat, df, type = 'response')
    m0 = predict(est_ctrl,  df, type = 'response')
  }
  # other predict methods go here
  if(ret_est == T){ # regular estimation
    if (est != 'rforest'){
      return(mean(m1 - m0))
    } else {
      # ranger returns other stuff too
      return(mean(m1[['predictions']] - m0[['predictions']]))
    }
  } else { # for aipw
    return(list(m0 = m0, m1 = m1))
  }
}


# %% ####################################################
#' Ingredients for AIPW estimator (naive)
#' @param meanfn name of outcome model
#' @param pscorefn name of pscore model
#' @param mean_fml formula for outcome model
#' @param psc_fml formula for pscore model
#' @param y outcome name
#' @param w treatment name
#' @param df dataframe
#' @export
fit_me = function(meanfn, pscorefn, mean_fml, psc_fml, y , w, df, ...){
  mean_fns = ate_reg(est = meanfn,   w = w, y = y, df = df, fml = mean_fml,
    ret_est = F, ...)
  pscores =      ate_ipw(est = pscorefn, w = w, y = y, df = df, fml = psc_fml,
    ret_est = F, ...)
  if (meanfn != 'rforest'){
    outdf = data.frame(df[[y]], df[[w]], mean_fns$m0, mean_fns$m1, pscores)
  } else {
    outdf = data.frame(df[[y]], df[[w]],
      mean_fns$m0$predictions, mean_fns$m1$predictions, pscores)
  }
  colnames(outdf) = c('y', 'w', 'm0', 'm1', 'ehat')
  outdf
}

# %% ####################################################
#' Ingredients for AIPW estimator (cross-fit)
#' @param meanfn name of outcome model
#' @param pscorefn name of pscore model
#' @param mean_fml formula for outcome model
#' @param psc_fml formula for pscore model
#' @param y outcome name
#' @param w treatment name
#' @param df dataframe
#' @param nfold number of folds (2 by default)
#' @export
fit_me_xf = function(meanfn, pscorefn, mean_fml, psc_fml, y , w, df, nfold = 2L, ...){
  nobs <- nrow(df) #number of observations
  foldid <- rep.int(1:nfold, times = ceiling(nobs/nfold))[sample.int(nobs)] # define folds indices
  I <- split(1:nobs, foldid) #split observation indices into folds
  # init m0, m1, ehat containers
  m0 <- m1  <- ehat <- rep(NA, nobs)
  #########################################
  # cross-fitting
  #########################################
  for(b in 1:length(I)){
    otherFolds = df[-I[[b]], ]; curFold = df[I[[b]], ]
    # treatment index
    treatInd = which(otherFolds[[w]] == 1)
    # compute model fits on -I fold
    m1mod = mhatter(estimator = meanfn,   fmla = mean_fml, y = y, df = otherFolds[treatInd,  ] )
    m0mod = mhatter(estimator = meanfn,   fmla = mean_fml, y = y, df = otherFolds[-(treatInd), ])
    emod  = ehatter(estimator = pscorefn, fmla = psc_fml,  w = w, df = otherFolds)
    # predict outcome in Ith fold
    if (meanfn %in% c('lasso', 'ridge')){
      m1[I[[b]]]     = predict(m1mod, model.matrix(mean_fml, curFold), s = 'lambda.min', type = 'response')
      m0[I[[b]]]     = predict(m0mod, model.matrix(mean_fml, curFold), s = 'lambda.min', type = 'response')
    } else if (meanfn == 'rforest'){
      m1[I[[b]]]     = predict(m1mod, model.matrix(mean_fml, curFold))[['predictions']]
      m0[I[[b]]]     = predict(m0mod, model.matrix(mean_fml, curFold))[['predictions']]
    } else if (meanfn == 'xgboost'){
      m1[I[[b]]]     = predict(m1mod, model.matrix(mean_fml, curFold))
      m0[I[[b]]]     = predict(m0mod, model.matrix(mean_fml, curFold))
    } else { # parametric models
      m1[I[[b]]]     = predict(m1mod, curFold, type = 'response')
      m0[I[[b]]]     = predict(m0mod, curFold, type = 'response')
    }
    # predict pscore in Ith fold
    if (pscorefn %in% c('lasso', 'ridge')){
      ehat[I[[b]]] = predict(emod, model.matrix(psc_fml, curFold), s = 'lambda.min', type = 'response')
    } else if (pscorefn == 'rforest'){
      ehat[I[[b]]] = predict(emod, model.matrix(psc_fml, curFold), type = 'response')[['predictions']]
    } else if (pscorefn == 'xgboost'){
      ehat[I[[b]]] = predict(emod, model.matrix(psc_fml, curFold))
    } else { # parametric models
      ehat[I[[b]]] = predict(emod, curFold, type = 'response')
    }
  }

  outdf = data.frame(df[[y]], df[[w]], m0, m1, ehat)
  colnames(outdf) = c('y', 'w', 'm0', 'm1', 'ehat')
  outdf
}

# %% ####################################################
#' compute AIPW
#' @param d output from fit_me function
#' @param psrange 2-vector of pscore to trim within
#' @return AIPW estimate of ATE
#' @export
ate_aipw = function(d, psrange = c(0, 1)){
  estsamp = d[d$ehat >= psrange[1] & d$ehat <= psrange[2],]
  est = with(estsamp,
    mean( (m1 + (w/ehat)*(y - m1) ) -
          (m0 + ( ((1 - w)/(1 - ehat)) * (y - m0)  ) )
      ))
  est
}
