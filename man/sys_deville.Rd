% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sys_deville.R
\name{sys_deville}
\alias{sys_deville}
\title{Deville's systematic}
\usage{
sys_deville(pik)
}
\arguments{
\item{pik}{A vector of inclusion probabilities.}
}
\value{
Return the selected indices in 1,2,...,N
}
\description{
This function implements a method to select a sample using the Deville's systmatic algorithm.
}
\examples{
set.seed(1)
pik <- c(0.2,0.5,0.3,0.4,0.9,0.8,0.5,0.4)
sys_deville(pik)
}
\references{
Deville, J.-C. (1998), Une nouvelle méthode de tirage à probabilité inégales. Technical Report 9804, Ensai, France.

 Chauvet, G. (2012), On a characterization of ordered pivotal sampling, Bernoulli, 18(4):1320-1340
}
