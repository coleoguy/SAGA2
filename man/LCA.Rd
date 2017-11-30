\name{LCA}
\alias{LCA}
\title{Model Averaged Line Cross Analysis}
\description{Analyze all possible genetic architecture models based on mean phenotypes from line cross data.}
\usage{LCA(data, SCS = "XY", parental = "calc", env = FALSE, model.sum = .95,
           max.models = 300000, drop.pars = NULL, max.pars = NULL, Cmatrix = NULL)}
\arguments{
  \item{data}{a data frame with the 7 columns:
  \itemize{
    \item cross: the name of each cohort e.g. P1, F1, BC1 etc.3
    \item mean: phenotype measure of the cohort
    \item SE: standard error of the cohort's mean phenotype
    \item sex: sex of the cohort \code{"M"}, \code{"F"}, \code{"E"}, and \code{"U"} representing male, female, equal ratio, unknown or unequal ratio are valid
    \item enviro: a numeric value indicating an environmental variable measured for each cohort - can be \code{NA} or \code{NULL}
    \item sire: the row of the user data that represents the sire of the current line
    \item dam: the row of the user data that represents the dam of the current line
    }
}
  \item{SCS}{A text string used to indicate the sex chromosome system \code{"XY"}, \code{"XO"}, \code{"ZW"}, \code{"ZO"}, or \code{"NSC"} (No Sex Chromosomes) are all valid.}
  \item{parental}{\code{"calc"} or \code{"obs"}
  \itemize{
     \item \code{"calc"}: a traditional calculated genetic maternal effect will be used
     \item \code{"obs"}: an observed phenotypic parental effect (maternal and paternal) will be used - see details for more information.}
     }
  \item{env}{\code{TRUE} or \code{FALSE} indicating whether to include E and GxE interactions.}
  \item{model.sum}{This is the sum of the probability of the models to be included}
  \item{max.models}{The maximum number of fitted models to return from the function.  This is included as an option to allow analysis of large model space on computers with limited RAM.}
  \item{drop.pars}{Vector of strings with the names of composite genetic effects to ignore.}
  }
  \item{max.pars}{Optional numeric value indicating the maximum number of parameters for equations evaluated.}
  \item{Cmatrix}{Optional matrix describing contribution of genetic effects in the cohorts being analyzed.}
  }
\details{Provides model averaged estimates of the contribution of composite genetic effects to the line means in line cross analysis experiments.  Using AICc models are given weights and these are used to construct a confidence model set that allows for parameter estimates and errors to include model selection uncertainty. (see Burnham and Anderson 2002).  The original paper contains a full discussion of the implementation (Blackmon and Demuth 2016).

In previous LCA work maternal effects have always been modeled as the autosomal additive or dominance status of the dam.  To the degree that the phenotype of interest is not controlled by either autosomal additive or dominance effects this approach will fail to capture the true nature of maternal effects.  As an alternative we offer the \code{parental} argument.  If supplied with the value \code{"obs"}  The phenotypes of all cohorts will be rescaled on an interval of -1 to 1. These values will then function as the expected contribution of either maternal or paternal effects in offspring from a cohort. 

In some experiments cohorts are combined to be used as parents in a subsequent cross. For instance an F1 and recipricol F1 might be combined and then backcrossed to each parent.  To indicate such a pooling of lines multiple rows of the user data can be indicated as either the sire or dam of a subsequent cohort.  Simple separate with a . and now spaces.  The example dataset shows what this should look like.

\strong{Example Data Set}

\tabular{lllllll}{
    \strong{cohort} \tab \strong{mean} \tab \strong{SE} \tab \strong{sex} \tab \strong{enviro} \tab \strong{sire} \tab \strong{dam}\cr
    P1 \tab 32 \tab .45 \tab E \tab 32 \tab 1 \tab 1\cr
    P2 \tab 16 \tab .65 \tab E \tab 32 \tab 2 \tab 2\cr
    F1 \tab 24 \tab 1.2 \tab E \tab 32 \tab 1 \tab 2\cr
    rF1 \tab 27 \tab 3.1 \tab E \tab 32 \tab 2 \tab 1\cr
    BC1a \tab 24 \tab 1.2 \tab E \tab 32 \tab 1 \tab 3.4\cr
    P1 \tab 29 \tab .01 \tab E \tab 29 \tab 6 \tab 6\cr
    P2 \tab 12 \tab .98 \tab E \tab 29 \tab 7 \tab 7\cr
    F1 \tab 22 \tab .8 \tab E \tab 29 \tab 6 \tab 7\cr
    rF1 \tab 21 \tab 2.4 \tab E \tab 29 \tab 7 \tab 6\cr
    rBC1a \tab 20 \tab .56 \tab E \tab 29 \tab 8.9 \tab 6
    }
 \itemize{
     \item Cohorts can be named in any way the user whishes but P1 and P2 are standard names that must be used for the two original strains.  If these are not present in the user data then the software will not be able to construct a C-matrix.
     \item Note that in the case of BC1a the dam is indicated to be a mix of both rF1 and F1 females.  The software will assume this is an equal mixture and is not able to analyze data with unequally pooled lines.
     \item This dataset also indicates that crosses were conducted in two environments (29 and 32). These will be rescaled appropraitely and treated as a potential environmental or gene by environmental effect on the phenotype being studied.
     }


}
\value{
Returns a "genarch" object which is a list with the following elements:
\item{all.models:}{a list containing the weighted least squares solution for 
all models tested.}
\item{best.models:}{a list containing the weighted least squares solution for 
all models in the confidence set.}
\item{best.eqns.w:}{a vector containing the model probabilities for each model in the confidence set.}
\item{estimates:}{a data frame containing Model Weighted Average for 
each parameter and its unconditional standard error.}
\item{daicc:}{a vector of the delta AICc scores for all models 
tested.}
\item{varimp:}{a data frame containing the variable importance scores for composite effects.}
\item{cmatrix:}{a data frame containing the contribution of possible genetic effects used in the analysis.}
}


\references{
Burnham, K. P., & Anderson, D. R. (2002). Model selection and multimodel 
inference: a practical information-theoretic approach. Springer.

Blackmon, H. and Demuth, J.P., (2016). An information‐theoretic approach to estimating the composite genetic effects contributing to variation among generation means: Moving beyond the joint‐scaling test for line cross analysis. Evolution, 70(2), pp.420-432.
}
\author{
Heath Blackmon and Jeffery P. Demuth}
\examples{
data(PH)
results <- LCA(PH)
}
