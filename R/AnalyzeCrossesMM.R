AnalyzeCrossesMM <- function(data, Cmatrix = "XY", 
                             env.factor=FALSE, even.sex = F,
                             max.models = 300000, model.sum = .95, max.pars = NULL){
  # lets check and make sure that people picked or supplied a cmatrix
  if(is.null(Cmatrix)) stop("Please supply or choose a Cmatrix")
  
  # if they are supplying the Cmatrix lets do a couple of basic checks
  if(!is.null(Cmatrix)){
    if(!is.vector(Cmatrix)){
      # is it even a matrix
      if(!is.matrix(Cmatrix)) stop("Your supplied c-matrix is not a matrix")
      # does it contain the cohorts we need
      if(!sum(data[,1] %in% Cmatrix[,1]) == nrow(data)){
        stop("The cohort IDs in your data don't match those in your c-matrix")
      }
    }
  }
  
  ## load the possible contributions to cohort means
  if(is.vector(Cmatrix)){
    if(Cmatrix == "XY" | Cmatrix == "XO" | Cmatrix == "X0"){
      Cmatrix <- read.csv(file = system.file("cmatrix.xy.csv", package = "SAGA2"), row.names=1)[, -1]
      Cmatrix <- Cmatrix[data[,1] , ]
      sex.dep <- c("sex", "Xa", "Xd", "Ya", "XaAa", "XaAd", "XdAa", "XdAd", 
                   "YaAa", "YaAd", "YaXa", "CaXa", "CaXd", "CaYa", "sex.Aa", 
                   "sex.Ad", "sex.Xa", "sex.Xd", "sex.Ca", "sex.Mea", "sex.Med",
                   "env.Xa", "env.Xd", "env.Ya")
      env.dep <- c("env", "env.Aa", "env.Ad", "env.Xa", "env.Xd", "env.Ya", "env.Ca", 
                   "env.Mea", "env.Med")
      y.dep   <- c("Ya", "YaAa", "YaAd", "YaXa", "CaYa", "sex.Ya", "env.Ya")
      if(even.sex == F){
        Cmatrix <- Cmatrix[, !colnames(Cmatrix) %in% sex.dep]
      }
      if(env.factor == F){
        Cmatrix <- Cmatrix[, !colnames(Cmatrix) %in% env.dep]
      }
      if(Cmatrix == "XO" | Cmatrix == "X0"){
        Cmatrix <- Cmatrix[, !colnames(Cmatrix) %in% y.dep]
      }
      
      
      
    } else if(Cmatrix == "ZW" | Cmatrix == "ZO" | Cmatrix == "Z0"){
      Cmatrix <- read.csv(file = system.file("cmatrix.zw.csv", package = "SAGA2"), row.names=1)[, -1]
      Cmatrix <- Cmatrix[data[,1] , ]
      sex.dep <- c("sex", "Za", "Zd", "Wa", "ZaAa", "ZaAd", "ZdAa", "ZdAd",
                   "WaAa", "WaAd", "ZaWa", "CaZa", "CaZd", "CaWa", "sex.Aa",
                   "sex.Ad", "sex.Za", "sex.Zd", "sex.Ca", "sex.Mea", "sex.Med",
                   "env.Za", "env.Zd", "env.Wa")
      env.dep <- c("env", "env.Aa", "env.Ad", "env.Za", "env.Zd", "env.Wa", "env.Ca", 
                   "env.Mea", "env.Med")
      w.dep   <- c("Wa", "WaAa", "WaAd", "ZaWa", "CaWa", "env.Wa")
      if(even.sex == F){
        Cmatrix <- Cmatrix[, !colnames(Cmatrix) %in% sex.dep]
      }
      if(env.factor == F){
        Cmatrix <- Cmatrix[, !colnames(Cmatrix) %in% env.dep]
      }
      if(Cmatrix == "ZO" | Cmatrix == "Z0"){
        Cmatrix <- Cmatrix[, !colnames(Cmatrix) %in% w.dep]
      }
      
      
    } else if(Cmatrix == "esd"){
      Cmatrix <- read.csv(file = system.file("cmatrix.esd.csv", package = "SAGA2"), row.names=1)[, -1]
      Cmatrix <- Cmatrix[data[,1] , ]
      sex.dep <- c("sex", "sex.Aa", "sex.Ad", "sex.Ca", "sex.Mea", "sex.Med")
      env.dep <- c("env", "env.Aa", "env.Ad", "env.Ca", "env.Mea", "env.Med")
      if(env.factor == T) Cmatrix[, 1] <- data[,2]
      # GxE
      Cmatrix[, 19] <- Cmatrix[, 1] * Cmatrix[, 4]
      Cmatrix[, 20] <- Cmatrix[, 1] * Cmatrix[, 5]
      Cmatrix[, 21] <- Cmatrix[, 1] * Cmatrix[, 6]
      Cmatrix[, 22] <- Cmatrix[, 1] * Cmatrix[, 7]
      Cmatrix[, 23] <- Cmatrix[, 1] * Cmatrix[, 8]
      if(even.sex == F){
        Cmatrix <- Cmatrix[, !colnames(Cmatrix) %in% sex.dep]
      }
      
      if(env.factor == F){
        Cmatrix <- Cmatrix[, !colnames(Cmatrix) %in% env.dep]
      }
      
    } else if(Cmatrix == "genomic"){
      print("Creating a genomic based cmatrix:")
      Cmatrix <- matrix(,nrow(data), 18)
      colnames(Cmatrix) <- c("env", "sex", "mean", "A", "D", "AA", "AD", "DD", 
                             "env.A", "env.D", "env.AA", "env.AD", "env.DD",
                             "sex.A", "sex.D", "sex.AA", "sex.AD", "sex.DD")
      if(env.factor == T) Cmatrix[, 1] <- data[,3]
      # TODO change to take vector of M/F
      if(sexed == T) Cmatrix[, 2] <- data[,4]
      # mean trait value
      Cmatrix[, 3] <- rep(1, times=nrow(Cmatrix))
      for(i in 1:nrow(Cmatrix)){
          s <- data[i, 1] # percent p1 genome
          h <- data[i, 2] # probability of being heterozygous p1/p2
          Cmatrix[i, 4] <- 1*s-(1*(1-s))
          Cmatrix[i, 5] <- h
          # digenic epistasis
          Cmatrix[i, 6] <- Cmatrix[i, 4] * Cmatrix[i, 4] 
          Cmatrix[i, 7] <- Cmatrix[i, 4] * Cmatrix[i, 5] 
          Cmatrix[i, 8] <- Cmatrix[i, 5] * Cmatrix[i, 5]
          # GxE
          Cmatrix[i, 9] <- Cmatrix[i, 1] * Cmatrix[i, 4]
          Cmatrix[i, 10] <- Cmatrix[i, 1] * Cmatrix[i, 5]
          Cmatrix[i, 11] <- Cmatrix[i, 1] * Cmatrix[i, 6]
          Cmatrix[i, 12] <- Cmatrix[i, 1] * Cmatrix[i, 7]
          Cmatrix[i, 13] <- Cmatrix[i, 1] * Cmatrix[i, 8]
          # GxSex
          Cmatrix[i, 14] <- Cmatrix[i, 2] * Cmatrix[i, 4]
          Cmatrix[i, 15] <- Cmatrix[i, 2] * Cmatrix[i, 5]
          Cmatrix[i, 16] <- Cmatrix[i, 2] * Cmatrix[i, 6]
          Cmatrix[i, 17] <- Cmatrix[i, 2] * Cmatrix[i, 7]
          Cmatrix[i, 18] <- Cmatrix[i, 2] * Cmatrix[i, 8]
        }
      sex.dep <- c("sex", "sex.A", "sex.D", "sex.AA", "sex.AD", "sex.DD")
      env.dep <- c("env", "env.A", "env.D", "env.AA", "env.AD", "env.DD")
      if(even.sex == F) Cmatrix[, !colnames(Cmatrix) %in% sex.dep]
      if(env.factor == F) Cmatrix[, !colnames(Cmatrix) %in% env.dep]
    } else {
      stop("Your selection for the c-matrix to use does not match any of the implemented options")
    }
  }
  
# lets remove variables that have no difference in lines
    cge.0 <- colnames(Cmatrix)[which(apply(Cmatrix, 2, var) == 0)]
    cge.0 <- cge.0[cge.0!="M"]
    red.Cmatrix <- Cmatrix
    if(length(cge.0) > 0){
      red.Cmatrix <- Cmatrix[, !colnames(Cmatrix) %in% cge.0]      
      print(paste(cge.0, "has no difference in expected line means and will not be considered in your model"))
    }

  
  # TODO we shouldn't be dropping higher order effects we should
  # let the user know that they are would receive equal support
  # and we should spread model uncertainty in accord with this
  #lets look for composite effects that are identical
  low.counter <- drop.counter <- vector()
  Aa <- which(colnames(red.Cmatrix) == "Aa")
  for(i in Aa:(ncol(red.Cmatrix)-1)){
    for(j in (i+1):ncol(red.Cmatrix)){
      if(sum(red.Cmatrix[,i] == red.Cmatrix[,j]) == nrow(red.Cmatrix)){
        drop.counter <- c(drop.counter, j)
        low.counter <-  c(low.counter, i)
      }
    }
  }
  #drop those composite effects that are equivelant of lower order simpler effects
  if(length(drop.counter)>0){
    leslie <- colnames(red.Cmatrix)[drop.counter]
    abbi <- colnames(red.Cmatrix)[low.counter]
    cat(paste("The following composite effects cannot be estimated with the line \n",
              "means available because they estimate identical quantities to \n",
              "lower order effects: \n",sep=""))
    for(i in 1:length(drop.counter)){
      cat(paste(leslie[i], "is being dropped because it equals", abbi, "\n", sep=""))
      red.Cmatrix <- red.Cmatrix[,-drop.counter]
    }
  }
  have.data <- paste(colnames(red.Cmatrix)[-1], collapse = ", ")
  cat(paste("The composite genetic effects that will be tested are: \n", 
              have.data, collapse = ", "), "\n\n")
  # calcualte the potential size of model space
  # the final -2 is because we will always be including the mean so we have
  # one less choice to make
  mod.space.size <- sum(choose(ncol(red.Cmatrix) -1 , 
                               1:(nrow(red.Cmatrix) - 2)))
  if(!is.null(max.pars)) mod.space.size <- sum(choose((ncol(red.Cmatrix) -1), 1:max.pars))
  # warn the user if the model space is very large
  if(mod.space.size > 5000){
    cat(paste("Since there are ", mod.space.size, " possible models this may take a bit:\n", sep=""))
  }
  # generate all possible models storing each matrix in a list
  pos.cols <- 2:ncol(red.Cmatrix)             # col that could be used
  eqns <- list()                              # store the eqns
  counter <- 1                                # index for eqns
  max.par <- nrow(red.Cmatrix) - 2            #
  # if a user has very many cohorts model space can become problematically large
  # however I think that actually very few datasets support >>large models with 
  # many important factors.  So one solution is simply to allow users to set a max
  # model size this makes things fairly easy to handle
  if(!is.null(max.pars)) max.par <- max.pars
  
  cat(paste("Generating Models"))
  if(length(pos.cols) < max.par){
    max.par <- length(pos.cols)
  }
  for(i in 1:max.par){                     # different number of par models
    cat(".")
    foo <- combn(pos.cols, i)              # all pos models with i variables
    # this loop just places the models generate with i variables into the list
    # of all possible models.  Models are described by the columns they include
    for(j in 1:ncol(foo)){
      eqns[[counter]] <- as.vector(foo[,j])
      counter <- counter + 1
    }
  }
  # just the setup for a small counter
  x <- 50
  if(length(eqns) > 1000) x <- 500
  if(length(eqns) > 10000) x <- 5000
  # now we test each model
  mod.results <- list()                  # stores glm results
  num.pars <- dev <- aic <- vector()     # stores various useful values
  # we need a counter because redundant models arrise.  These originate because
  # some components will have high covariance depending on the lines included
  # in the dataset.  The glm function automatically throws these variables
  # resulting in fitting the same model more than once.
  counter <- 1
  for(i in 1:length(eqns)){
    # generate the matrix for the current model
    test.mat <- as.matrix(red.Cmatrix[, c(1, eqns[[i]])])
    # fit the model weight is equal to the inverse of the square of the SE
    temp.mod <- glm(data[, 1] ~ test.mat, weights = data[, 2] ^ - 2)
    # this if statement will bypass a model with a singularity
    # 1 NA will be generated for the line mean any additional are sign of sing.
    if(sum(is.na(temp.mod$coef)) < 2){
      # name model results as eqns
      mod.results[[counter]] <- temp.mod
      names(mod.results)[counter] <- i
      # record the number of parameters in the model
      num.pars[counter] <- length(mod.results[[counter]]$coefficients) - 1
      # record the residual deviances
      dev[counter] <- mod.results[[counter]]$dev
      # record the AIC of the models
      aic[counter] <- mod.results[[counter]]$aic
      counter <- counter + 1
    }
    if(i %% x == 0) cat(paste("\n", i))
  }
  ## need to report the number of models thrown out due to 
  ## high covariance ~ singularity
  if(i > counter){
  cat(paste("\n", i - (counter - 1), 
            " models were removed due to high covariances \n",
              "or linear relationships between predictor variables.  \n", "The remaining ", 
              counter - 1, " models have been evaluated.\n\n", sep = ""))
  }
  # in the unrealistic situation where there was a model that predicted the
  # data perfectly we would get -Inf for the AIC should only be an issue in 
  # simulated data... for these purposes lets just plug in something that is 
  # equal to the lowest AIC value for models in that same parameter range
  for(i in 1:length(aic)){
    if(aic[i] == -Inf){
      aic[i] <- sort(unique(aic))[2] -1
      print("One of the models fit the data perfectly... Hopefully you are using
            simulated data.  If not please investigate becuase this shouldn't be
            happening.")
    }
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
  new.vars <- matrix(0,(ncol(red.Cmatrix)-1),2)
  new.vars[,1] <- colnames(red.Cmatrix)[2:ncol(red.Cmatrix)]
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
  cat(paste("\nAICc weights were used to select the minimum number of models ",
              "whose weights sum \nto greater than ", 
            model.sum * 100, "% this model set includes ", length(best.models), 
              " model(s)\n", sep = ""))
   #lets calculate variable importance
   #which equations are being used
   best.eqns <- eqns[as.numeric(names(best.models))]
   best.eqns.w <- waic[sort(best.models.ind)]
   #lets tell the user the models being included
   #for(i in 1:length(best.eqns)){
   #     cat(paste(colnames(red.Cmatrix)[c(1, best.eqns[[i]])], collapse = ", "), 
   #      "  waic = ",good.model.waics[i],"\n")
   #  }
  # now we need to print the model weighted averages and SE
  # lets make a matrix of the calculated values under each model
  par.est <- matrix(0, length(best.eqns), ncol(red.Cmatrix) + 2)
  colnames(par.est) <- c('eqn', colnames(red.Cmatrix), 'mw')
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
  results <- matrix(, 2, ncol(red.Cmatrix))
  colnames(results) <- colnames(red.Cmatrix)
  row.names(results) <- c('Model Weighted Average', 
                          'Unconditional Standard Error')
  results[1, ] <- par.est[nrow(par.est), 2:(ncol(par.est) - 1)]
  results[2, ] <- var.est[nrow(var.est), 2:(ncol(var.est) - 1)]
  
  ## prepare the results to be returned to the user
  final.results <- list()
  mod.names <- list()
  foo2 <- colnames(red.Cmatrix)
  mod.ind <- as.numeric(names(mod.results))
  for(i in 1:length(mod.results)){
    mod.names[[i]] <- paste(foo2[eqns[[mod.ind[i]]]], sep="", collapse=", ")
  }
  names(mod.results) <- mod.names
  if(length(mod.results) > max.models){
    mod.results <- mod.results[order(waic, decreasing = T)[1:max.models]]
    daicc <- daicc[order(daicc, decreasing =F)[1:max.models]]
  }
  final.results[[1]] <- mod.results
  final.results[[2]] <- results
  final.results[[3]] <- daicc
  final.results[[4]] <- new.vars
  names(final.results) <- c("models", "estimates", "daicc", "varimp")
  class(final.results) <- "genarch"
  return(final.results)
  }
  