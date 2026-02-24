#' Model Averaged Line Cross Analysis
#'
#' Analyze all possible genetic architecture models based on mean phenotypes
#' from line cross data. Uses an information-theoretic approach with AICc-based
#' model selection to produce model averaged estimates of composite genetic
#' effects (CGEs) and unconditional standard errors.
#'
#' @param data A data frame with 7 columns:
#'   \describe{
#'     \item{cross}{Character. Name of each cohort (e.g., "P1", "P2", "F1",
#'       "rF1", "BC1a"). P1 and P2 must be present for the two parental
#'       strains.}
#'     \item{mean}{Numeric. Mean phenotype measure of the cohort.}
#'     \item{SE}{Numeric. Standard error of the cohort mean.}
#'     \item{sex}{Character. Sex of the cohort: `"M"` (male), `"F"` (female),
#'       `"E"` (equal ratio), or `"U"` (unknown/unequal).}
#'     \item{enviro}{Numeric. Environmental variable for each cohort. Can be
#'       `NA` if not applicable.}
#'     \item{sire}{Integer. Row number of the sire in the data frame.}
#'     \item{dam}{Integer. Row number of the dam in the data frame.}
#'   }
#' @param SCS Character. Sex chromosome system: `"XY"`, `"XO"`, `"ZW"`,
#'   `"ZO"`, or `"NSC"` (no sex chromosomes). Default is `"XY"`.
#' @param parental Character. Method for calculating parental effects:
#'   `"calc"` uses traditional calculated genetic maternal effects;
#'   `"obs"` uses observed phenotypic parental effects rescaled to \[-1, 1\].
#'   Default is `"calc"`.
#' @param env Logical. Whether to include environmental and gene-by-environment
#'   interactions. Default is `FALSE`.
#' @param model.sum Numeric between 0 and 1. Cumulative probability threshold
#'   for the confidence set of models. Default is `0.95`.
#' @param max.models Integer. Maximum number of fitted model objects to store
#'   and return. Only affects memory usage; all models are still evaluated
#'   internally. Default is `300000`.
#' @param drop.pars Character vector. Names of CGEs to exclude from the
#'   analysis (e.g., `c("Ca", "Mea")`). Default is `NULL`.
#' @param keep.pars Character vector. Names of CGEs to include in the analysis.
#'   If specified, only these CGEs (plus the mean) will be tested. Default is
#'   `NULL`.
#' @param max.pars Integer. Maximum number of parameters allowed in any single
#'   model equation. Useful for reducing computation time with large model
#'   spaces. Default is `NULL` (determined automatically).
#' @param Cmatrix Optional matrix. A user-supplied C-matrix describing the
#'   expected contribution of CGEs to each cohort. If `NULL`, the C-matrix
#'   is constructed automatically from the breeding design. Default is `NULL`.
#' @param ret.all Logical. Whether to return fitted model objects for all
#'   evaluated models. Setting to `TRUE` increases memory usage substantially.
#'   Default is `FALSE`.
#' @param messages Logical. Whether to print progress messages during the
#'   analysis. Default is `TRUE`.
#' @param Mepi Logical. Whether to include interactions between maternal
#'   effects and standard genetic effects. Experimental. Default is `FALSE`.
#'
#' @details
#' The analysis proceeds by:
#' 1. Constructing a C-matrix of expected CGE contributions (or using a
#'    user-supplied one).
#' 2. Fitting all possible weighted least squares models.
#' 3. Calculating AICc and Akaike weights for each model.
#' 4. Constructing a confidence set of models.
#' 5. Computing model-averaged parameter estimates and unconditional
#'    standard errors following Burnham and Anderson (2002).
#'
#' When `parental = "obs"`, cohort phenotypes are rescaled to \[-1, 1\] and
#' used as the expected contribution of maternal or paternal effects, which
#' can better capture parental effects not mediated through autosomal additive
#' or dominance mechanisms.
#'
#' @return A `"genarch"` object (a list) with the following elements:
#' \describe{
#'   \item{all.models}{A list of fitted `glm` model objects for all evaluated
#'     models (or `NULL` if `ret.all = FALSE`).}
#'   \item{best.models}{A list of fitted `glm` model objects for models in
#'     the confidence set.}
#'   \item{best.eqns.w}{Numeric vector of Akaike weights for models in the
#'     confidence set.}
#'   \item{estimates}{A 2-row matrix with model-weighted average parameter
#'     estimates (row 1) and unconditional standard errors (row 2) for each
#'     CGE.}
#'   \item{daicc}{Numeric vector of delta AICc scores for all evaluated
#'     models.}
#'   \item{varimp}{A 2-column matrix with CGE names (column 1) and variable
#'     importance scores (column 2).}
#'   \item{cmatrix}{The C-matrix used in the analysis, with cross names
#'     prepended.}
#' }
#'
#' @references
#' Blackmon, H. and Demuth, J.P. (2016). An information-theoretic approach to
#' estimating the composite genetic effects contributing to variation among
#' generation means: Moving beyond the joint-scaling test for line cross
#' analysis. *Evolution*, 70(2), 420--432.
#'
#' Burnham, K.P. and Anderson, D.R. (2002). *Model selection and multimodel
#' inference: a practical information-theoretic approach*. Springer.
#'
#' @seealso [plot.genarch()] for plotting results, [VisModelSpace()] for
#'   visualizing model space, [EvaluateModel()] for single model evaluation,
#'   [plotObserved()] for traditional LCA plots.
#'
#' @export
#' @examples
#' data(BanInf)
#' result <- LCA(BanInf, SCS = "NSC", messages = FALSE)
#' result$estimates
#' plot(result)
LCA <- function(data,
                SCS = "XY",
                parental = "calc",
                env = FALSE,
                model.sum = .95,
                max.models = 300000,
                drop.pars = NULL,
                keep.pars = NULL,
                max.pars = NULL,
                Cmatrix = NULL,    # user supplied cmatrix
                ret.all = FALSE,   # return solution for all models
                messages = TRUE,
                Mepi = FALSE){     # maternal effect epistatic interactions

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
  data$cross <- toupper(data$cross)
  ### if no custom matrix is supplied build a cmatrix based
  ### on the user data and arguments supplied
  if(is.null(Cmatrix)){
    Cmatrix <- PrepareCmatrix(user.data = data,
                              SCS = SCS, env = env,
                              drop.pars = drop.pars,
                              keep.pars = keep.pars,
                              parental = parental,
                              Mepi = Mepi)
  }

  # remove CGEs in the cmatrix that can't be analyzed
  cmat.temp <- CleanCmatrix(Cmatrix, messages=messages)
  Cmatrix <- cmat.temp[[1]]
  crosses <- cmat.temp[[2]]

  ### report the composite genetic effects being explored
  have.data <- paste(colnames(Cmatrix)[-1], collapse = ", ")
  if(messages == T) cat(paste("The composite genetic effects that will be tested are: \n",
                              have.data, "\n", collapse = ", "))

  ### calculate the potential size of model space
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
