\name{AICtoMW}
\alias{AICtoMW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
converts AIC or AICc to model weights}
\description{
This function generates model weights based on either delta AIC scores or delta AICc scores.}
\usage{
AICtoMW(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector containing delta AIC or delta AICc scores}
  }
\details{
This function returns a numeric vector containing model weights

}
\author{
Heath Blackmon and Jeffery P. Demuth}
\examples{
# load data
data(PH)
# analyze it
results <- AnalyzeCrossesMM(PH)
# calculate model weights
AICtoMW(results[[3]])
}
