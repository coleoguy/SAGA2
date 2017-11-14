\name{AnalyzeCrossesMM}
\alias{AnalyzeCrossesMM}
\title{Multimodel Analysis of Line Cross Data}
\description{Analyze all possible genetic architecture models based on mean phenotypes from line cross data.}
\usage{AnalyzeCrossesMM(data, Cmatrix = "XY", model.sum = .95,
                 max.models = 300000, even.sex = F, graph=F, 
                 cex.axis=1, cex.names=1, cex.main=1, max.pars = NULL)}
\arguments{
  \item{data}{a data frame with the first three columns:
1) id of the cohort this must mach the coefficient row of the c-matrix
2) mean phenotype measure of the cohort
3) Standard error of the cohort's mean phenotype
}
  \item{Cmatrix}{A text string used to select the mid-parent scaled c-matrix to be used in the analysis.  Included options are "XY", "XO", "ZW", "ZO", or "esd".  If you need a different matrix you can supply it here as well.}
  \item{model.sum}{This is the sum of the probability of the models to be included}
  \item{max.models}{The maximum number of fitted models to return from the function.  This is included as an option to allow analysis of large model space on computers with limited RAM.}
  \item{even.sex}{A logical by default it is false. It should be set as true if either sexed cohorts are included or if mixed sex cohorts are included but have equal numbers of males and females.}
  \item{max.pars}{Optional parameter limiting the size of the equations evaluated.}
  \item{graph}{Logical indicating whether a plot of results should be produced}
  \item{cex.axis}{expansion factor for numeric axis labels.}
  \item{cex.names}{expansion factor for CGE labels.}
  \item{cex.main}{expansion factor for main title.}
  }
\details{Provides model averaged estimates of the contribution of composite genetic effects to the line means in line cross analysis experiments.  Using AICc models are given weights and these are used to construct a confidence model set that allows for parameter estimates and errors to include model selection uncertainty. (see Burnham and Anderson 2002).  The vignette contains a lengthy discussion of 
}
\value{
Returns a "genarch" object which is a list with the following elements:
\item{models:}{a list containing the weighted least squares solution for 
all models tested.}
\item{estimates:}{a data frame containing Model Weighted Average for 
each parameter and its unconditional standard error.}
\item{daicc:}{a vector of the delta AICc scores for all models 
tested.}
\item{varimp:}{a data frame containing the variable importance scores for composite effects.}
}
\references{
Burnham, K. P., & Anderson, D. R. (2002). Model selection and multimodel 
inference: a practical information-theoretic approach. Springer.
}
\author{
Heath Blackmon and Jeffery P. Demuth}
\examples{
data(PH)
results <- AnalyzeCrossesMM(PH)
}
\seealso{
\code{\link{VisModelSpace}}
}
