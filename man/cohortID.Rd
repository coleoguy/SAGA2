\name{cohortID}
\alias{cohortID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Display a lookup table for cohort IDs
}
\usage{
cohortID(sexed=T)
}
\arguments{
  \item{sexed}{if True then will return all available cohorts if False then returns only mixed sex cohorts}
}
\description{
A data frame allowing the user to lookup ID numbers for each cohort.  Parents are listed in parentheses and sire x dam.
}

\author{
Heath Blackmon and Jeffery P. Demuth}
\examples{
IDs <- cohortID()
}
