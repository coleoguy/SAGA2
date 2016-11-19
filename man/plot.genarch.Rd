\name{plot}
\alias{plot.genarch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
provides a custom plot of a genarch object}
\description{
Plot a publication quality figure of a genarch object with many built-in customization such as only a subset with variable importances over a chosen cut-off.  If a cut-off is chosen that would make it so only a single bar would be plotted the min.vi argument will be automatically reduced to include at least the two CGEs with the highest variable importance.}
\usage{
\method{plot}{genarch}(x, min.vi = .5, main=NULL, 
             cex.axis=1, cex.names=1, cex.main=1, 
             maxval=NULL, minval=NULL, viridis=F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{This is the observed data used in the function AnalyzeCrossesMM}
  \item{min.vi}{This is a "genarch" object returned by the function AnalyzeCrossesMM}
  \item{main}{title for the Y axis}
  \item{cex.axis}{expansion factor for Y axis}
  \item{cex.names}{expansion factor for CGE names on x axis}
  \item{cex.main}{expansion factor for main title}
  \item{maxval}{custom maximum value for a plot}
  \item{minval}{custom minimum value for a plot}
  \item{viridis}{logical indicating if the viridis color palette should be used}
  \item{\dots}{further arguments passed to plot}
}
\author{
Heath Blackmon and Jeffery P. Demuth}
\examples{
data(PH)
results <- AnalyzeCrossesMM(PH)
plot(results)
}
