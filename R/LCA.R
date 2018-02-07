LCA <- function(data,
                SCS = "XY",
                parental = "calc",
                env = FALSE,
                model.sum = .95,
                max.models = 300000,
                drop.pars = NULL,
                max.pars = NULL,
                Cmatrix = NULL,
                ret.all = FALSE,
                messages = TRUE){
  
  ### lets deal with cross names being treated as factors
  if(is.factor(data$cross)){
    data$cross <- unlist(lapply(data$cross, as.character))
  }
  
  ### lets deal with F in sex being treated as FALSE
  if(is.logical(data$sex)){
    data$sex <- rep("F", length(data$sex))
  }
  
  # validate the incoming arguments and data
  validateData(SCS = SCS, user.data = data, Cmatrix = Cmatrix, messages = messages)
  
  ### if no custom matrix is supplied build a cmatrix based
  ### on the user data and arguments supplied
  if(is.null(Cmatrix)){
    Cmatrix <- PrepareCmatrix(user.data = data,
                              SCS = SCS, env = env,
                              drop.pars = drop.pars,
                              parental = parental)
  }
  
  # remove CGEs in the cmatrix that can't be analyzed
  cmat.temp <- CleanCmatrix(Cmatrix, messages=messages)
  Cmatrix <- cmat.temp[[1]]
  crosses <- cmat.temp[[2]]
  
  ### report the composite genetic effects being explored
  have.data <- paste(colnames(Cmatrix)[-1], collapse = ", ")
  if(messages == T) cat(paste("The composite genetic effects that will be tested are: \n",
                              have.data, "\n", collapse = ", "))
  
  ### calcualte the potential size of model space
  mod.space.size <- sum(choose(ncol(Cmatrix) -1 , 1:(nrow(Cmatrix) - 2)))
  if(!is.null(max.pars)) mod.space.size <- sum(choose((ncol(Cmatrix) -1), 1:max.pars))
  if(messages == T){
    if(mod.space.size > 5000){
      cat(paste("Since there are", mod.space.size,
                "possible models this may take a bit:\n"))
    }
  }
  
  ### analyze the data based on the cmatrix
  result <- AnalyzeModels(data = data,
                          Cmatrix = Cmatrix,
                          crosses = crosses,
                          SCS = SCS,
                          parental = parental,
                          env = env,
                          model.sum = model.sum,
                          max.models = max.models,
                          max.pars = max.pars,
                          ret.all = ret.all,
                          messages = messages)
  return(result)
}
