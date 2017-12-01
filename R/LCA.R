LCA <- function(data,
                SCS = "XY",
                parental = "calc",
                env = FALSE,
                model.sum = .95,
                max.models = 300000,
                drop.pars = NULL,
                max.pars = NULL,
                Cmatrix = NULL){
  # lets deal with cross names being treated as factors
  if(is.factor(data$cross)){
    data$cross <- unlist(lapply(data$cross, as.character))
  }
  
  # lets deal with F in sex being treated as FALSE
  if(is.logical(data$sex)){
    data$sex <- rep("F", length(data$sex))
  }
  
  # validate the incoming arguments and data
  validateData(SCS, user.data = data, Cmatrix)
  # if no custom matrix is supplied build a cmatrix based
  # on the user data and arguments supplied
  if(is.null(Cmatrix)){
    Cmatrix <- PrepareCmatrix(user.data = data,
                              SCS = SCS, env = env,
                              drop.pars = drop.pars,
                              parental = parental)
  }
  # remove CGEs in the cmatrix that can't be analyzed
  cmat.temp <- CleanCmatrix(Cmatrix)
  Cmatrix <- cmat.temp[[1]]
  crosses <- cmat.temp[[2]]
  # analyze the data based on the cmatrix
  result <- AnalyzeModels(data = data,
                          Cmatrix = Cmatrix,
                          crosses = crosses,
                           SCS = SCS,
                           parental = parental,
                           env = env,
                           model.sum = model.sum,
                           max.models = max.models,
                           max.pars = max.pars)
  return(result)
}
