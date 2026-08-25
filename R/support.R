# These are internal functions used in SAGA 2.0

######### rescales 0,1 to -1,1 #########
ReScale <- function(x) {
  2 * x - 1
}
########################################

######### rescale vector to -1,1 #######
ReScale2 <- function(x) {
  min.val <- min(x)
  max.val <- max(x)
  z <- ((2 * x - 2 * min.val) / (max.val - min.val))-1
  return(z)
}
########################################

######### validate data ################
# Validates the user supplied data and arguments. Because R passes
# arguments by value, any fixes applied here (SCS synonyms, missing
# environmental values, SE column renaming) must be returned to the
# caller; validateData returns a list with the corrected data and SCS.
validateData <- function(SCS, user.data, Cmatrix, messages){
  # if they are supplying the Cmatrix lets do a couple of basic checks
  if(!is.null(Cmatrix)){
    if(!is.vector(Cmatrix)){
      # is it even a matrix
      if(!is.matrix(Cmatrix)) stop("Your supplied c-matrix is not a matrix")
      # does it contain the cohorts we need
      if(!sum(user.data[,1] %in% Cmatrix[,1]) == nrow(user.data)){
        stop("The cross names in your data don't
             match those in your user supplied c-matrix")
      }
      # the first CGE column must be the all-1s mean column; the rest of
      # the package assumes this column is present in every model
      if(!isTRUE(all(as.numeric(Cmatrix[, 2]) == 1))){
        stop("The second column of a user supplied c-matrix must be the
             mean column and contain only 1s.")
      }
    }
  }
  # validate and fix problems with SCS supplied by user
  if (SCS == "X0") SCS <- "XO"
  if (SCS == "Z0") SCS <- "ZO"
  if (!SCS %in% c("XY", "XO", "ZW", "ZO", "NSC")) {
    stop(
      "The sex chromosome system supplied should be either
    XY, XO, ZW, ZO, or NSC.  NSC (no sex chromosome system)
    should be used if you are unsure of the sex chromosome
    system present."
    )
  }
  # Here is a validation test that the sex state is either
  # male, female, equal, or unknown/unequal
  if (!all(user.data$sex %in% c("M", "F", "E", "U"))) {
    stop(
      "Sex of each line should be either M, F, E, or
    U indicating either males, females, equal ratio,
    or unequal"
    )
  }
  # Test that phenotypes are numeric values
  if (!is.numeric(user.data$mean)) {
    stop("Phenotypes must be numeric values")
  }
  # Test that standard errors are numeric values
  colnames(user.data)[3] <- toupper(colnames(user.data)[3])
  if (!is.numeric(user.data$SE)) {
    stop("Standard Errors must be numeric values")
  }
  if (any(is.na(user.data$SE))) {
    stop(paste("Standard Errors contain NA values in row(s):",
               paste(which(is.na(user.data$SE)), collapse = ", ")))
  }
  # If the user supplies no environmental values fill with 1
  if(all(is.na(user.data$enviro))) user.data$enviro <- 1

  # Test that environment is between -1 and 1
  if (!is.numeric(user.data$enviro)){
    stop("Environment values need to be numeric values.  If you have discrete
       environmental states use -1 and 1.")
  }
  # Test that every sire and dam value is a valid row number
  if(is.null(Cmatrix)){
    n <- NROW(user.data)
    bad <- which(is.na(user.data$sire) | is.na(user.data$dam) |
                   user.data$sire < 1 | user.data$dam < 1 |
                   user.data$sire > n | user.data$dam > n |
                   user.data$sire != round(user.data$sire) |
                   user.data$dam != round(user.data$dam))
    if (length(bad) > 0) {
      stop(paste("Sire and dam values must be integer row numbers between 1
    and the number of rows in your data (", n, "). Problem row(s): ",
                 paste(bad, collapse = ", "), sep = ""))
    }
  }
  # return the corrected data and SCS to the caller
  return(list(data = user.data, SCS = SCS))
}


######### End of validation testing ####

######### Clean the Cmatrix #############
CleanCmatrix <- function(Cmatrix, messages=FALSE){
  # lets pull out crosses
  crosses <- Cmatrix[, 1]
  Cmatrix <- Cmatrix[, -1]
  # lets remove variables that have no difference in lines
  Cmatrix <- Cmatrix[, c(1, which(apply(Cmatrix, 2, var) != 0))]
  #lets look for composite effects that are identical
  # (duplicated() on the transposed matrix flags every column that is
  # identical to an earlier column and safely handles minimal c-matrices)
  drop.counter <- which(duplicated(t(Cmatrix)))
  drop.counter <- drop.counter[drop.counter >= 2]
  #drop those composite effects that are equivelant of lower order simpler effects
  if(length(drop.counter)>0){
    badCGEs <- paste(colnames(Cmatrix)[drop.counter], sep=", ", collapse=", ")
    if(messages==TRUE){
      cat(paste("The following composite effects cannot be estimated with the line \n",
                "means available because they estimate identical quantities to \n",
                "lower order effects: \n", badCGEs, "\n\n", sep=""))
    }
    Cmatrix <- Cmatrix[,-drop.counter]
  }
  return(list(Cmatrix, crosses))
}
######### End of Cmatrix cleaning ######
