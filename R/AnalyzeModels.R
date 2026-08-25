AnalyzeModels <- function(data, 
                          SCS, 
                          Cmatrix, 
                          crosses, 
                          parental, 
                          env, 
                          model.sum, 
                          max.models, 
                          max.pars,
                          ret.all = TRUE, 
                          messages){

  ##### generate all possible models storing each matrix in a list
  # col that could be used
  pos.cols <- 2:ncol(Cmatrix)
  # store the eqns
  eqns <- list() 
  # index for eqns
  counter <- 1
  # the maximum allowed parameters
  max.par <- nrow(Cmatrix) - 2           
  # allow user to overide max.par  
  if(!is.null(max.pars)) max.par <- max.pars
  
  if(messages == TRUE) cat(paste("Generating Models"))
  if(length(pos.cols) < max.par){
    max.par <- length(pos.cols)
  }
  
  # different number of par models
  for(i in 1:max.par){                     
    if(messages == TRUE) cat(".")
    foo <- combn(pos.cols, i)              
    # Models are described by a vector of the columns they include
    for(j in 1:ncol(foo)){
      eqns[[counter]] <- as.vector(foo[,j])
      counter <- counter + 1
    }
  }
  
  # just the setup for a small counter
  if(length(eqns) <= 1000) x <- 50
  if(length(eqns) > 1000) x <- 500
  if(length(eqns) > 10000) x <- 5000
  
  # if the user supplies the matrix make sure it is numeric
  if(!is.null(Cmatrix)){
    if(is.character(Cmatrix[1, 1])){
      class(Cmatrix) <- "numeric"
    }
  }
  # convert the c-matrix to a numeric matrix once rather than once per model
  mm <- as.matrix(Cmatrix)
  storage.mode(mm) <- "numeric"

  # weights are constant across models so compute them once;
  # weight is equal to the inverse of the square of the SE
  SEs <- data$SE
  if(any(SEs == 0)){
    SEs[SEs == 0] <- .000000000001
  }
  weights <- SEs ^ -2
  n.obs <- nrow(data)
  sum.log.w <- sum(log(weights))

  # now we test each model.  For speed the exhaustive scoring pass uses
  # lm.wfit (identical weighted least squares fits to glm but without the
  # per-model formula/model-frame overhead); glm objects are refit below
  # only for the models that are returned to the user.
  n.eqns <- length(eqns)
  keep <- logical(n.eqns)
  num.pars <- dev <- aic <- rep(NA_real_, n.eqns)

  for(i in 1:n.eqns){
    # generate the matrix for the current model
    X <- mm[, c(1, eqns[[i]]), drop = FALSE]
    fit <- lm.wfit(x = X, y = data$mean, w = weights)
    # skip models with a singularity (rank deficient fits); these arise
    # because some components have high covariance
    if(fit$rank == ncol(X)){
      keep[i] <- TRUE
      # record the number of parameters in the model
      num.pars[i] <- fit$rank
      # record the residual deviance
      dv <- sum(weights * fit$residuals ^ 2)
      dev[i] <- dv
      # record the AIC of the model (matches stats::glm for
      # gaussian family fits: family aic + 2 * rank)
      aic[i] <- n.obs * (log(2 * pi * dv / n.obs) + 1) + 2 -
        sum.log.w + 2 * fit$rank
    }
    if(messages == TRUE) if(i / x == round(i / x)) cat(paste("\n", i))
  }

  # keep only the models that could be fit
  kept <- which(keep)
  counter <- length(kept)
  num.pars <- num.pars[kept]
  dev <- dev[kept]
  aic <- aic[kept]

  ## need to report the number of models thrown out due to 
  ## high covariance ~ singularity
  if(messages == TRUE){
    if(n.eqns > counter){
      cat(paste("\n", n.eqns - counter, 
                " models were removed due to high covariances \n",
                "or linear relationships between predictor variables.  \n", "The remaining ", 
                counter, " models have been evaluated.\n\n", sep = ""))
    }
  }
  
  # second order AIC (AICc); K counts all estimated parameters
  # including the residual variance (sigma^2)
  K <- num.pars + 1
  aicc <- ifelse(n.obs - K - 1 > 0,
                 aic + (2 * K * (K + 1)) / (n.obs - K - 1),
                 Inf)
  daicc <- aicc - min(aicc)
  
  # this code correctly produces akaike weights 
  waic <- (exp(-.5 * daicc)) / (sum(exp(-.5 * daicc)))
  
  # calculate variable importance: sum the weights of every model that
  # contains each CGE
  cge.names <- colnames(mm)[2:ncol(mm)]
  memb <- unlist(eqns[kept])
  w.rep <- rep(waic, lengths(eqns[kept]))
  imp <- rowsum(w.rep, memb)
  vi <- rep(0, length(cge.names))
  vi[as.numeric(rownames(imp)) - 1] <- imp[, 1]
  new.vars <- matrix(vi, ncol = 1,
                     dimnames = list(cge.names, "Importance"))
  
  # lets calculate the 95% probability set of models
  o <- order(waic, decreasing = TRUE, na.last = NA)
  n.conf <- which(cumsum(waic[o]) >= model.sum)[1]
  if(is.na(n.conf)) n.conf <- length(o)
  best.models.ind <- o[1:n.conf]
  if(messages == TRUE){
    cat(paste("\nAICc weights were used to select the minimum number of models ",
              "whose weights sum \nto greater than ", 
              model.sum * 100, "% this model set includes ", length(best.models.ind), 
              " model(s)\n", sep = ""))
  }
  
  # refit the returned models with glm so users receive standard fitted
  # model objects; this is restricted to the confidence set (and, when
  # requested, ret.all models) to keep the exhaustive pass fast
  fit.glm <- function(eqn.ind){
    test.mat <- mm[, c(1, eqns[[eqn.ind]]), drop = FALSE]
    glm(data$mean ~ test.mat, weights = weights)
  }
  
  # which equations are being used and their weights; best.eqns.w is in
  # the same (decreasing weight) order as best.models
  best.eqn.ind <- kept[best.models.ind]
  best.models <- lapply(best.eqn.ind, fit.glm)
  names(best.models) <- best.eqn.ind
  best.eqns <- eqns[best.eqn.ind]
  best.eqns.w <- waic[best.models.ind]
  
  # now we need the model weighted averages and unconditional SEs;
  # all bookkeeping is numeric with membership tracked in a logical matrix
  p <- ncol(mm)
  nb <- length(best.models)
  memb.mat <- matrix(FALSE, nb, p, dimnames = list(names(best.models),
                                                   colnames(mm)))
  est.mat <- var.mat <- matrix(0, nb, p,
                               dimnames = list(names(best.models),
                                               colnames(mm)))
  for(i in seq_len(nb)){
    cols <- c(1, best.eqns[[i]])
    memb.mat[i, cols] <- TRUE
    # drop the aliased coefficient by its NA status rather than position
    cf <- coef(best.models[[i]])
    est.mat[i, cols] <- cf[!is.na(cf)]
    var.mat[i, cols] <- diag(vcov(best.models[[i]], complete = FALSE))
  }
  # recalculate model waic to sum to 1
  mw <- best.eqns.w / sum(best.eqns.w)
  
  # calculate the model weighted parameter estimates
  avg <- colSums(est.mat * mw)
  
  # now lets calculate the unconditional standard errors as proposed in
  # burnham and anderson 2002 pg 162
  uncond.se <- rep(0, p)
  for(i in seq_len(p)){
    if(avg[i] != 0){
      uncond.se[i] <- sum(mw * sqrt(var.mat[, i] +
                                      (est.mat[, i] - avg[i]) ^ 2))
    }
  }
  
  # now lets make a table with the stuff we want
  results <- matrix(NA_real_, 2, p)
  colnames(results) <- colnames(mm)
  row.names(results) <- c('Model Weighted Average', 
                          'Unconditional Standard Error')
  results[1, ] <- avg
  results[2, ] <- uncond.se
  
  ## prepare the results to be returned to the user
  foo2 <- colnames(mm)
  keep.for.all <- seq_len(counter)
  if(counter > max.models){
    keep.for.all <- order(waic, decreasing = TRUE)[1:max.models]
    daicc <- daicc[order(daicc, decreasing = FALSE)[1:max.models]]
  }
  mod.results <- NULL
  if(ret.all == TRUE){
    all.eqn.ind <- kept[keep.for.all]
    mod.results <- lapply(all.eqn.ind, fit.glm)
    names(mod.results) <- vapply(all.eqn.ind, function(ix){
      paste(foo2[eqns[[ix]]], sep = "", collapse = ", ")
    }, character(1))
  }
  final.results <- list(all.models = mod.results,
                        best.models = best.models,
                        best.eqns.w = best.eqns.w,
                        estimates = results,
                        daicc = daicc,
                        varimp = new.vars,
                        cmatrix = cbind(crosses, Cmatrix))
  class(final.results) <- "genarch"
  return(final.results)
}