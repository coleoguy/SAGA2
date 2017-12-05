PrepareCmatrix <- function(user.data, 
                           SCS, 
                           parental, 
                           drop.pars, 
                           env, 
                           messages=T) {
  
  ##### Scale environmental factors #####
  if(!all(min(user.data$enviro) != -1, max(user.data$enviro) != 1)){
    # if there is only one value we dont need to do anything
    if(length(unique(user.data$enviro)) > 1){
      if(messages == T){
        print("Scaling environmental variable to be on a scale of -1, 1")
      }
      user.data$enviro <- ReScale2(user.data$enviro)
    }
  }
  ##### Environmental factors are now scaled #####
  
  
  ##### Now we create an empty dataframe for the P-matrix #####
  pmatrix <- as.data.frame(matrix(, nrow(user.data), 7))
  # the f and m in the col names below indicate the value
  # in a female and a male respectively
  switch(SCS,
         "XY" = , "XO" = , "NSC" = colnames(pmatrix) <-
           # cross    %p1 autosome  %fem p1X   male %p1X  Y source  Mito source   phenotype measure
           c("cross", "p1a",        "fp1x",    "mp1x",    "y",      "c",          "pheno"),
         "ZW" = , "ZO" = colnames(pmatrix) <-
           c("cross", "p1a", "fp1z", "mp1z", "w", "c", "pheno")
  )
  ##### empty P-matrix done#####
  
  ##### Now we create an empty dataframe for the C-matrix #####
  switch(SCS,
         "XY" = cmatrix <- as.data.frame(matrix(, nrow(user.data), 11)),
         "ZW" = cmatrix <- as.data.frame(matrix(, nrow(user.data), 11)),
         "XO" = cmatrix <- as.data.frame(matrix(, nrow(user.data), 10)),
         "ZO" = cmatrix <- as.data.frame(matrix(, nrow(user.data), 10)),
         "NSC" = cmatrix <- as.data.frame(matrix(, nrow(user.data), 8))
  )
  
  # The col names below have the same meanings as in the orginal SAGA
  # paper and application.
  switch(SCS,
         "XY" = colnames(cmatrix)  <- c("cross", "M", "Aa", "Ad",
                                        "Xa", "Xd", "Ya",
                                        "Ca", "Mea", "Med", "Env"),
         "ZW" = colnames(cmatrix)  <- c("cross", "M", "Aa", "Ad",
                                        "Za", "Zd", "Wa",
                                        "Ca", "Mea", "Med", "Env"),
         "XO" = colnames(cmatrix)  <- c("cross", "M", "Aa", "Ad", "Xa", "Xd",
                                        "Ca", "Mea", "Med", "Env"),
         "ZO" = colnames(cmatrix)  <- c("cross", "M", "Aa", "Ad", "Za", "Zd",
                                        "Ca", "Mea", "Med", "Env"),
         "NSC" = colnames(cmatrix) <- c("cross", "M", "Aa", "Ad",
                                        "Ca", "Mea", "Med", "Env")
  )
  
  # the col names below (Meo and Peo) are new additions for SAGA2
  # and indicate observed phenotypic scores rescaled on a -1 to 1 range
  if (parental == "obs") {
    colnames(cmatrix)[colnames(cmatrix) %in% c("Mea", "Med")] <-
      c("Meo", "Peo")
  }
  cmatrix$cross <- user.data$cross
  ##### Lets check for pooled parents ######
  #first check to see if there are any
  # pool.count <- 0
  # if(sum(user.data$sire %% 1) > 0 |
  #    sum(user.data$dam %% 1) > 0){
  #   parents <- c(user.data$sire, user.data$dam)
  #   pool.types <- unique(parents[parents %% 1 != 0])
  #   pool.count <- length(pool.types)
  # }
  # if(pool.count > 0){
  #   # I've added this bit to deal with cases where we have more than one P1 or P2
  #   # for instance we may have males and females seperate or we may have common
  #   # garden experiments where species are measured in different environments
  #   # create the pool name
  #   parents <- as.numeric(strsplit(as.character(pool.types), split=".", fixed=T)[[1]])
  #   pname <- paste(user.data$cross[parents[1]], user.data$cross[parents[2]], sep=".",collapse="")
  # }
  #get cross names
  pmatrix$cross <- user.data$cross
  # first we get the P1 and P2 rows
  prefill <- which(pmatrix$cross %in% c("P1","P2"))
  # now we will fill these rows
  for(i in prefill){
    switch(user.data$cross[i],
           "P1" = pmatrix[i, 2:7] <- c(1, 1, 1,  1,  1, 
                                       user.data$mean[i]),
           "P2" = pmatrix[i, 2:7] <- c(0, 0, 0, -1, -1, 
                                       user.data$mean[i]))
  }
  # get unfilled columns
  to.fill <- which(!complete.cases(pmatrix))
  for(i in to.fill){
    # this chunk is for rows that dont have any pooled parents
      pmatrix$p1a[i]   <- (pmatrix$p1a[user.data$sire[i]] + 
                             pmatrix$p1a[user.data$dam[i]]) / 2.0
      pmatrix$c[i]     <- pmatrix$c[user.data$dam[i]]
      pmatrix$pheno[i] <- user.data$mean[i]
      if(SCS == "XY" || SCS == "XO" || SCS == "NSC"){
        pmatrix$fp1x[i]  <- (pmatrix$mp1x[user.data$sire[i]] + 
                               pmatrix$fp1x[user.data$dam[i]]) / 2.0
        pmatrix$mp1x[i]  <- pmatrix$fp1x[user.data$dam[i]]
        pmatrix$y[i]     <- pmatrix$y[user.data$sire[i]]
      }
      if(SCS == "ZW" || SCS == "ZO"){
        pmatrix$mp1z[i]  <- (pmatrix$mp1z[user.data$sire[i]] + 
                               pmatrix$fp1z[user.data$dam[i]]) / 2.0
        pmatrix$fp1z[i]  <- pmatrix$mp1z[user.data$sire[i]]
        pmatrix$w[i]     <- pmatrix$w[user.data$dam[i]]
      }
  }
  ##### pmatrix filled#####
  
  
  
  
  ##### fnx to fill cmatrix #####
  FillCmat <- function(cmatrix, pmat, user.data) {
    # to minimize the number of lines with really
    # long code we will user litte composite
    # effect calculator functions
    #Aa
    pc <- function(pmatrix, user.data, i) {
      ReScale((pmatrix$p1a[user.data$sire[i]] +
                 pmatrix$p1a[user.data$dam[i]]) / 2)
    }
    #Ad
    hc <- function(pmatrix, user.data, i) {
      (pmatrix$p1a[user.data$sire[i]] +
         pmatrix$p1a[user.data$dam[i]]) -
        2 * pmatrix$p1a[user.data$sire[i]] *
        pmatrix$p1a[user.data$dam[i]]
    }
    # Xa in males
    pxmc <- function(pmatrix, user.data, i) {
      ReScale(pmatrix$fp1x[user.data$dam[i]])
    }
    # Xa in females
    pxfc <- function(pmatrix, user.data, i) {
      ReScale((pmatrix$fp1x[user.data$dam[i]] +
                 pmatrix$mp1x[user.data$sire[i]]) / 2)
    }
    # Xd in females
    hxfc <- function(pmatrix, user.data, i) {
      (pmatrix$fp1x[user.data$dam[i]] +
         pmatrix$mp1x[user.data$sire[i]]) -
        2 * pmatrix$fp1x[user.data$dam[i]] *
        pmatrix$mp1x[user.data$sire[i]]
    }
    # Za in males
    pzmc <- function(pmatrix, user.data, i) {
      ReScale((pmatrix$fp1z[user.data$dam[i]] +
                 pmatrix$mp1z[user.data$sire[i]]) / 2)
    }
    # Za in females
    pzfc <- function(pmatrix, user.data, i) {
      ReScale(pmatrix$mp1z[user.data$dam[i]])
    }
    # Zd in males
    hzmc <- function(pmatrix, user.data, i) {
      (pmatrix$fp1z[user.data$dam[i]] +
         pmatrix$mp1z[user.data$sire[i]]) -
        2 * pmatrix$fp1z[user.data$dam[i]] *
        pmatrix$mp1z[user.data$sire[i]]
    }
    #now we want to loop through all our cols and rows
    for (i in 1:nrow(cmatrix)) {
      for (j in 1:ncol(cmatrix)) {
        switch(colnames(cmatrix)[j],
               "M"   = cmatrix[i, j] <- 1,
               "Aa"  = cmatrix[i, j] <- pc(pmatrix, user.data, i),
               "Ad"  = cmatrix[i, j] <- hc(pmatrix, user.data, i),
               "Mea" = cmatrix[i, j] <- cmatrix$Aa[user.data$dam[i]],
               "Med" = cmatrix[i, j] <- cmatrix$Ad[user.data$dam[i]],
               "Meo" = cmatrix[i, j] <- pmatrix$pheno[user.data$dam[i]],
               "Peo" = cmatrix[i, j] <- pmatrix$pheno[user.data$sire[i]],
               "Env" = cmatrix[i, j] <- user.data$enviro[i],
               "Ca"  = cmatrix[i, j] <- pmatrix$c[i]
        )
        if (user.data$sex[i] == "F") {
          switch(colnames(cmatrix)[j],
                 "Xa" = cmatrix[i, j] <- pxfc(pmatrix, user.data, i),
                 "Xd" = cmatrix[i, j] <- hxfc(pmatrix, user.data, i),
                 "Za" = cmatrix[i, j] <- pzfc(pmatrix, user.data, i),
                 "Wa" = cmatrix[i, j] <- pmatrix$w[user.data$dam[i]]
          )
        }
        if (user.data$sex[i] == "M") {
          switch(colnames(cmatrix)[j],
                 "Xa" = cmatrix[i, j] <- pxmc(pmatrix, user.data, i),
                 "Za" = cmatrix[i, j] <- pzmc(pmatrix, user.data, i),
                 "Zd" = cmatrix[i, j] <- hzmc(pmatrix, user.data, i),
                 "Ya" = cmatrix[i, j] <- pmatrix$y[i]
          )
        }
        if (user.data$sex[i] == "E") {
          switch(colnames(cmatrix)[j],
                 "Xa" = cmatrix[i, j] <-
                   (pxfc(pmatrix, user.data, i) + pxmc(pmatrix, user.data, i)) / 2,
                 "Xd" = cmatrix[i, j] <- hxfc(pmatrix, user.data, i) / 2,
                 "Za" = cmatrix[i, j] <-
                   (pzfc(pmatrix, user.data, i) + pzmc(pmatrix, user.data, i)) / 2,
                 "Zd" = cmatrix[i, j] <- hzmc(pmatrix, user.data, i) / 2,
                 "Ya" = cmatrix[i, j] <- pmatrix$y[user.data$sire[i]] / 2,
                 "Wa" = cmatrix[i, j] <- pmatrix$w[user.data$dam[i]] / 2
          )
        }
      }
    }
    drops <- names(cmatrix)[apply(cmatrix,FUN=function(x) all(is.na(x)), MARGIN = 2)]
    cmatrix <- cmatrix[, !names(cmatrix) %in% drops]
    return(cmatrix)
  }
  ##### end of fnx to fill cmatrix #####
  
  ##### fill the basic cmatrix #####
  cmatrix <- FillCmat(cmatrix, pmat, user.data)
  ##### fill the basic cmatrix #####
  
  ##### add env to drop if there are no diffs #####
  if(length(unique(user.data$enviro)) < 2){
    if(is.null(drop.pars)) drop.pars <- "Env"
    if(!is.null(drop.pars)){
      if(!"Env" %in% drop.pars){
        drop.pars <- c(drop.pars, "Env")
      }
    }
  }
  ##### drop user specified CGEs #####
  
  ##### drop user specified CGEs #####
  if(!is.null(drop.pars)){
    x <- which(colnames(cmatrix) %in% drop.pars)
    cmatrix <- cmatrix[, -x]
  }
  
  if(env == FALSE){
    cmatrix <- cmatrix[, names(cmatrix) != "Env"]
  }  
  
  
  
  ##### end of user restriction of CGEs #####
  
  ##### expand the cmatrix to include epistatic #####
  v <- 3:ncol(cmatrix)
  v <- which(!colnames(cmatrix)[3:ncol(cmatrix)] %in% c("Mea", "Med", "Meo", "Peo")) + 2
  n <- length(v)
  new.cols <- combinations(n = n, r = 2, v = v, repeats.allowed = T)
  for(i in 1:nrow(new.cols)){
    x <- new.cols[i, 1]
    y <- new.cols[i, 2]
    temp.col <- cmatrix[,x] * cmatrix[,y]
    cmatrix <- cbind(cmatrix, temp.col)
    names(cmatrix)[ncol(cmatrix)] <- paste(names(cmatrix)[x], names(cmatrix)[y], sep="")
  }
  # drop EnvEnv that makes no sense
  cmatrix <- cmatrix[, names(cmatrix) != "EnvEnv"]
  ##### end of cmatrix expansion #####
  
  ##### return the gold to the user #####
    return(cmatrix)
    return(pmatrix)
  ##### we should be done #####
}
