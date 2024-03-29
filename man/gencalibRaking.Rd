% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gencalibRaking}
\alias{gencalibRaking}
\title{Generalized calibration using raking ratio}
\usage{
gencalibRaking(Xs, Zs, d, total, q, max_iter = 500L, tol = 1e-09)
}
\arguments{
\item{Xs}{A matrix of calibration variables.}

\item{Zs}{A matrix of instrumental variables with same dimension as Xs.}

\item{d}{A vector, the initial weights.}

\item{total}{A vector that represents the initial weights.}

\item{q}{A vector of positive value that account for heteroscedasticity.}

\item{max_iter}{An integer, the maximum number of iterations. Default = 500.}

\item{tol}{A scalar that represents the tolerance value for the algorithm. Default = 1e-9.}
}
\value{
A vector, the value of the g-weights.
}
\description{
This function is inspired by the function \code{\link[sampling:calib]{calib}} of the package sampling. It computes the g-weights of the calibration estimator.
}
\details{
More details on the different calibration methods can be read in Tillé Y. (2020).
}
\references{
Tillé, Y. (2020). \emph{Sampling and estimation from finite populations}. Wiley, New York
}
\author{
Raphaël Jauslin \email{raphael.jauslin@unine.ch}
}
