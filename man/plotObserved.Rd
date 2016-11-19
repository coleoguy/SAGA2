\name{plotObserved}
\alias{plotObserved}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
provides the traditional line cross analysis plot.}
\description{
plot the observed results from a line cross analysis experiment as well as expectation from a simple additive model.}
\usage{
plotObserved(data, results, pch, col, xlab, ylab, main, SE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{This is the observed data used in the function AnalyzeCrossesMM}
  \item{results}{This is a "genarch" object returned by the function AnalyzeCrossesMM}
  \item{pch}{optional - a vector to specify the symbol to be used for observed line means if the vector is shorter than the number of lines it will be recycled.}
  \item{col}{optional - a vector to specify the color of symbol to be used for observed line means if the vector is shorter than the number of lines it will be recycled.}
  \item{xlab}{optional - title for the X axis}
  \item{ylab}{optional - title for the Y axis}
  \item{main}{optional - main title for the graph}
  \item{SE}{include standard error bars; default is TRUE}
}
\author{
Heath Blackmon and Jeffery P. Demuth}
\examples{
data(PH)
plotObserved(data=PH)
}
