\name{EvaluateModel}
\alias{EvaluateModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate Single Model}
\description{
This function allows the user to explore a specific model in a "genarch" object}
\usage{
EvaluateModel(data, model, cex.axis=1, cex.names=1, cex.main=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
This is a "genarch" object returned by the function AnalyzeCrossesMM}
  \item{model}{
This is the number corresponding to the model the user wishes to investigate.}
  \item{cex.axis}{
expansion factor for numeric axis labels.}
  \item{cex.names}{
expansion factor for name labels.}
  \item{cex.main}{
expansion factor for main title.}
}
\author{
Heath Blackmon and Jeffery P. Demuth}
\examples{
data(PH)
results <- AnalyzeCrossesMM(PH)
EvaluateModel(results, model=11)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
