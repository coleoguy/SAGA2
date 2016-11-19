\name{data}
\alias{data}
\alias{ban.osa}
\alias{per.inf}
\alias{dar.bho}
\alias{sin.cro}
\alias{PH}
\alias{SL}
\alias{SR}

\docType{data}
\title{empirical datasets}
\description{

ban.osa: Productivity of crosses involving Tribolium castaneum from Ecuador and Japan [1]

dar.bho: Productivity of crosses involving Tribolium castaneum from Tanzania and India [1]

per.inf: Productivity of crosses involving Tribolium castaneum from Peru and Portugal [1]

sin.cro: Productivity of crosses involving Tribolium castaneum from Malaysia and Croatia [1]

PH: Plant height data for Nicotiana rustica [2]

SL: Sperm length data for male Drosophila mojavensis line cross  analysis [3]

SR: Sperm receptacle length data for female Drosophila mojavensis line cross analysis [3]
}
\usage{data(ban.osa)}
\value{A data frame with rows for each cohort and 3 columns:
  \item{id}{ a numeric vector }
  \item{mean}{ a numeric vector }
  \item{SE}{ a numeric vector }
  }
\source{
[1] Demuth, J.P., Evolution of Hybrid Incompatibility in the beetle Tribolium Castaneum, in Biology. 2004, Indiana University: Bloomington. p. 152.

[2] Lynch, M., and B. Walsh, 1998 Lynch, M., & Walsh, B. (1998). Genetics and analysis of quantitative traits. Sinauer Associates, Inc., Sunderland, Massachusetts.

[3] Miller, G. T., Starmer, W. T., & S. Pitnick 2003. Quantitative genetic analysis of among-population variation in sperm and female sperm-storage organ length in Drosophila mojavensis. Genetical research, 81(03), 213-220.

}
\examples{
data(ban.osa)
data(dar.bho)
data(per.inf)
data(sin.cro)
data(PH)
data(SL)
data(SR)
}
\keyword{datasets}