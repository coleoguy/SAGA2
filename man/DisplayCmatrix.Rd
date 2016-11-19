\name{DisplayCmatrix}
\alias{DisplayCmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Display a matrix of composite genetic effects
}
\usage{
DisplayCmatrix(table = "XY")
}
\arguments{
  \item{table}{a character string specifying the sex determination system.  Current options are "XY", ZW", "esd"}
}
\description{
A data frame representing the c-matrix scaled to the midparent. Column names describe the composite effects.  The argument "esd" should be used for species that do not have sex chromosomes.}
\value{
  \item{M}{ mean }
  \item{Aa}{ autosomal additive }
  \item{Ad}{ autosomal dominance }
  \item{Xa}{ X chromosome additive }
  \item{Xd}{ X chromosome dominance }
  \item{Ya}{ Y chromosome additive }
  \item{Ca}{ cytotype additive }
  \item{Mea}{ maternal effect additive }
  \item{Med}{ maternal effect dominance }
  \item{AaAa}{ autosomal additive by additive epistasis }
  \item{AaAd}{ autosomal additive by dominance epistasis }
  \item{AdAd}{ autosomal dominance by dominance epistasis }
  \item{XaAa}{ X chromosome additive by autosomal additive epistasis }
  \item{XaAd}{ X chromosome additive by autosomal dominance epistasis }
  \item{XdAa}{ X chromosome dominance by autosomal additive epistasis }
  \item{XdAd}{ X chromosome dominance by autosomal dominance epistasis }
  \item{YaAa}{ Y chromosome additive by autosomal additive epistasis }
  \item{YaAd}{ Y chromosome additive by autosomal dominance epistasis }
  \item{YaXa}{ Y chromosome additive by X additive epistasis }
  \item{CaAa}{ cytotype additive by autosomal additive epistasis }
  \item{CaAd}{ cytotype additive by autosomal dominance epistasis }
  \item{CaXa}{ cytotype additive by X chromosome additive epistasis }
  \item{CaXd}{ cytotype additive by X chromosome dominance epistasis }
  \item{CaYa}{ cytotype additive by Y chromosome additive epistasis }
}

\author{
Heath Blackmon and Jeffery P. Demuth}
\examples{
cmat <- DisplayCmatrix("XY")
}
