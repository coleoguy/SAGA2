AnalyzeModels <- function(data, 
                          SCS, 
                          Cmatrix, 
                          crosses, 
                          parental, 
                          env, 
                          model.sum, 
                          max.models, 
                          max.pars,
                          ret.all, 
                          messages){
  
  # report the composite genetic effects being explored
  have.data <- paste(colnames(Cmatrix)[-1], collapse = ", ")
  if(messages == T){
    cat(paste("The composite genetic effects that will be tested are: \n", 
            have.data, collapse = ", "), "\n\n")
  }
  
  # calcualte the potential size of model space
  # the final -2 is because we will always be including the mean so we have
  # one less choice to make
  mod.space.size <- sum(choose(ncol(Cmatrix) -1 , 
                               1:(nrow(Cmatrix) - 2)))
  if(!is.null(max.pars)) mod.space.size <- sum(choose((ncol(Cmatrix) -1), 1:max.pars))
  # warn the user if the model space is very large
  if(messages == T){
    if(mod.space.size > 5000){
      cat(paste("Since there are ", mod.space.size, " possible models this may take a bit:\n", sep=""))
    }
  }
  
  # generate all possible models storing each matrix in a list
  pos.cols <- 2:ncol(Cmatrix)             # col that could be used
  eqns <- list()                              # store the eqns
  counter <- 1                                # index for eqns
  max.par <- nrow(Cmatrix) - 2            #
  
  # if a user has very many cohorts model space can become problematically large
  # however I think that actually very few datasets support >>large models with 
  # many important factors.  So one solution is simply to allow users to set a max
  # model size this makes things fairly easy to handle
  if(!is.null(max.pars)) max.par <- max.pars
  if(messages == T) cat(paste("Generating Models"))
  if(length(pos.cols) < max.par){
    max.par <- length(pos.cols)
  }
  # different number of par models
  for(i in 1:max.par){                     
    if(messages == T) cat(".")
    foo <- combn(pos.cols, i)              
    # all pos models with i variables
    # this loop just places the models generate with i variables into the list
    # of all possible models.  Models are described by the columns they include
    for(j in 1:ncol(foo)){
      eqns[[counter]] <- as.vector(foo[,j])
      counter <- counter + 1
    }
  }
  # just the setup for a small counter
  if(length(eqns) < 101) x <- 50
  if(length(eqns) > 100) x <- 50
  if(length(eqns) > 1000) x <- 500
  if(length(eqns) > 10000) x <- 5000
  # now we test each model

  # We need to preallocate these variables
  # mod.results, num.pars, dev, aic
  mod.results <- vector(mode = "list", length = length(eqns))
  num.pars <- dev <- aic <- vector(length=length(eqns))
  
  # we need a counter because redundant models arrise.  These originate because
  # some components will have high covariance depending on the lines included
  # in the dataset.  The glm function automatically throws these variables
  # resulting in fitting the same model more than once.
  counter <- 0
  for(i in 1:length(eqns)){
    # generate the matrix for the current model
    test.mat <- as.matrix(Cmatrix[, c(1, eqns[[i]])])
    # fit the model weight is equal to the inverse of the square of the SE
    temp.mod <- glm(data$mean ~ test.mat, weights = data$SE ^ - 2)
    # this if statement will bypass a model with a singularity
    # 1 NA will be generated for the line mean any additional are sign of sing.
    if(sum(is.na(temp.mod$coef)) < 2){
      counter <- counter + 1
      # name model results as eqns
      mod.results[[counter]] <- temp.mod
      names(mod.results)[counter] <- i
      # record the number of parameters in the model
      num.pars[counter] <- length(mod.results[[counter]]$coefficients) - 1
      # record the residual deviances
      dev[counter] <- mod.results[[counter]]$dev
      # record the AIC of the models
      aic[counter] <- mod.results[[counter]]$aic
    }
    if(messages == T) if(i / x == round(i / x)) cat(paste("\n", i))
  }
  
  # Get rid of excess preallocation
  mod.results <- mod.results[1:counter]
  num.pars <- num.pars[1:counter]
  dev <-  dev[1:counter]
  aic <- aic[1:counter]
  
  
  
  
  ## need to report the number of models thrown out due to 
  ## high covariance ~ singularity
  if(messages == T){
    if(i > counter){
      cat(paste("\n", i - (counter - 1), 
                " models were removed due to high covariances \n",
                "or linear relationships between predictor variables.  \n", "The remaining ", 
                counter - 1, " models have been evaluated.\n\n", sep = ""))
    }
  }
  # in the unrealistic situation where there was a model that predicted the
  # data perfectly we would get -Inf for the AIC should only be an issue in 
  # simulated data... for these purposes lets just plug in something that is 
  # equal to the lowest AIC value for models in that same parameter range
  if(min(aic) == -Inf){
    aic[aic == -Inf] <- sort(unique(aic))[2] -1
  }
  # calculate aicc and delta aicc
  aicc <- aic + (((2 * num.pars) * (num.pars + 1)) / 
                   (nrow(data) - num.pars))
  daicc <- aicc - min(aicc)
  # this code correctly produces akaike weights 
  waic <- (exp(-.5 * daicc)) / (sum(exp(-.5 * daicc)))
  new.waic <- waic
  new.waic.names <- eqns[as.numeric(names(mod.results))]
  # so now we have a copy of the waics to play with
  new.vars <- matrix(0,(ncol(Cmatrix)-1),2)
  new.vars[,1] <- colnames(Cmatrix)[2:ncol(Cmatrix)]
  for(i in 1:nrow(new.vars)){
    for(j in 1:length(new.waic)){
      if((i+1) %in% new.waic.names[[j]]){
        new.vars[i,2] <- as.numeric(new.vars[i,2]) + new.waic[j]
      }
    }
  }
  # lets calculate the 95% probability set of models
  best.models <- list()
  counter <- i <- 0
  good.model.waics <- vector()
  while(counter < model.sum){
    i <- i + 1
    counter <- counter + waic[order(waic, decreasing = T, na.last = F)][i]
    good.model.waics[i] <- waic[order(waic, decreasing = T, na.last = F)][i]
  }
  best.models.ind <- order(waic, decreasing = T, na.last = F)[1:i]
  best.models <- mod.results[best.models.ind]
  if(messages == T){
    cat(paste("\nAICc weights were used to select the minimum number of models ",
              "whose weights sum \nto greater than ", 
              model.sum * 100, "% this model set includes ", length(best.models), 
              " model(s)\n", sep = ""))
  }
  
  #lets calculate variable importance
  #which equations are being used
  best.eqns <- eqns[as.numeric(names(best.models))]
  best.eqns.w <- waic[sort(best.models.ind)]
  
  # now we need to print the model weighted averages and SE
  # lets make a matrix of the calculated values under each model
  par.est <- matrix(0, length(best.eqns), ncol(Cmatrix) + 2)
  colnames(par.est) <- c('eqn', colnames(Cmatrix), 'mw')
  par.est[, 1] <- names(best.models)
  par.est[, 2] <- 1
  
  # now we need a 1 or 0  if the parameter is in the eqn
  for(i in 1:nrow(par.est)){
    bar <- as.numeric(par.est[i, 1])
    par.est[i, eqns[[bar]] + 1] <- 1
  }
  
  # now replace 1's with the parameter estimate for each variable
  for(i in 1:nrow(par.est)){
    bar <- best.models[[i]]$coefficients[-2]
    counter <- 0
    for(j in 2:ncol(par.est)){
      if(par.est[i, j] == 1){
        counter <- counter + 1
        par.est[i, j] <- bar[counter]
      }
    }
  }
  
  # add in the aicw
  # best.models has eqn lookup in waic
  names(waic) <- names(mod.results)
  for(i in 1:nrow(par.est)){
    par.est[i, ncol(par.est)] <- waic[names(waic) == par.est[i, 1]]
  }
  
  # recalculate model waic to sum to 1
  par.est[, 'mw'] <- as.numeric(par.est[, 'mw']) / 
    sum(as.numeric(par.est[, 'mw']))
  
  # calculate the model weighted parameter estimates
  par.est <- rbind(par.est, rep(0,ncol(par.est)))
  for(i in 2:(ncol(par.est)-1)){
    par.est[nrow(par.est), i] <- sum(as.numeric(par.est[1:nrow(par.est)-1, i]) * 
                                       as.numeric(par.est[1:nrow(par.est)-1, 'mw']))
  }
  par.est[nrow(par.est), 1] <- 'mw.avg'
  # now lets calculate the unconditional variances as proposed in burnham and
  # anderson 2002 pg 162
  # so lets duplicate table par.est to use to fill in our values
  var.est <- par.est
  # lets loop through models first with i ... the rows
  for(i in 1:(nrow(var.est) - 1)){
    counter <- 0
    mod.vars <- diag(vcov(best.models[[which(names(best.models) == 
                                             var.est[i, 1])]]))
    # now lets loop through parameters with j ... the columns
    for(j in 2:(ncol(var.est) - 1)){
      if(var.est[i, j] != 0){
        counter <- counter + 1
        var.est[i,j] <- mod.vars[counter] 
      }
    }
  }
  
  # we now have all of the required variables for the uncond. variance est.
  for(i in 2:(ncol(var.est) - 1)){
    if(as.numeric(par.est[nrow(par.est), i]) != 0){
      foo <- vector()
      for(j in 1:(nrow(var.est) - 1)){
        foo[j] <- as.numeric(par.est[j, 'mw']) * 
          sqrt(as.numeric(var.est[j, i]) + 
                 ((as.numeric(par.est[j, i]) - 
                     as.numeric(par.est[nrow(par.est), i])) ^ 2))
      }
      var.est[nrow(var.est), i] <- sqrt(sum(foo) ^ 2)
    }
  }
  
  # now lets make a table with the stuff we want
  results <- matrix(, 2, ncol(Cmatrix))
  colnames(results) <- colnames(Cmatrix)
  row.names(results) <- c('Model Weighted Average', 
                          'Unconditional Standard Error')
  results[1, ] <- par.est[nrow(par.est), 2:(ncol(par.est) - 1)]
  results[2, ] <- var.est[nrow(var.est), 2:(ncol(var.est) - 1)]
  
  ## prepare the results to be returned to the user
  final.results <- list()
  mod.names <- list()
  foo2 <- colnames(Cmatrix)
  mod.ind <- as.numeric(names(mod.results))
  for(i in 1:length(mod.results)){
    mod.names[[i]] <- paste(foo2[eqns[[mod.ind[i]]]], sep="", collapse=", ")
  }
  names(mod.results) <- mod.names
  if(length(mod.results) > max.models){
    mod.results <- mod.results[order(waic, decreasing = T)[1:max.models]]
    daicc <- daicc[order(daicc, decreasing =F)[1:max.models]]
  }
  if(ret.all == TRUE) final.results[[1]] <- mod.results
  if(ret.all == FALSE) final.results[[1]] <- NULL
  final.results[[2]] <- best.models
  final.results[[3]] <- best.eqns.w
  final.results[[4]] <- results
  final.results[[5]] <- daicc
  final.results[[6]] <- new.vars
  final.results[[7]] <- cbind(crosses, Cmatrix)
  names(final.results) <- c("all.models", "best.models", 
                            "best.eqns.w", "estimates", 
                            "daicc", "varimp",
                            "cmatrix")
  class(final.results) <- "genarch"
  return(final.results)
}