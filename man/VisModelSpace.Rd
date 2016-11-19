\name{VisModelSpace}
\alias{VisModelSpace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize Model Space}
\description{
This function plots squares for every model evaluated.  The color of the square
is a function of the AICc weight of that model.  Model numbers are printed 
on top of each square allowing the user to investigate any models of interest
using the function EvaluateModel}
\usage{
VisModelSpace(data, cex.u = 3, cex.mtext = 1, cex.leg = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{This is a "genarch" object returned by the function AnalyzeCrossesMM}
  \item{cex.u}{expansion factor for the boxes representing the models}
  \item{cex.mtext}{expansion factor for the x and y axis labels}
  \item{cex.leg}{expansion factor for the legend}
}
\author{
Heath Blackmon and Jeffery P. Demuth}
\examples{
data(PH)
results <- AnalyzeCrossesMM(PH)
VisModelSpace(results, cex.u=3)
}
