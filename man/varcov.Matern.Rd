% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{varcov.Matern}
\alias{varcov.Matern}
\title{Matern covariance function}
\usage{
varcov.Matern(X, nu = 5/2)
}
\arguments{
\item{X}{vector to compute the covariance matrix}

\item{nu}{scalar; smoothness parameter of the Matern covariance matrix. Default at 5/2.}
}
\value{
a matrix of the dimension n x n where n is the length of X.
}
\description{
This function computes the Matern covariance function for a vector X.
}
\keyword{CondCop}
