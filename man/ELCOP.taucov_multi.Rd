% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{ELCOP.taucov_multi}
\alias{ELCOP.taucov_multi}
\title{Nonparametric Bayesian estimate of the conditional Kendall's tau through the empirical likelihood.}
\usage{
ELCOP.taucov_multi(
  uu,
  x1,
  x2,
  nsim = 10000,
  bdw,
  kern = "gauss",
  n1grid = 30,
  n2grid = 30
)
}
\arguments{
\item{uu}{n x 2 matrix of pseudo-observations. n is the sample size.}

\item{x1}{vector of n values of covariate 1}

\item{x2}{vector of n values of covariate 2}

\item{nsim}{number of Monte Carlo simulations. Default at 10,000.}

\item{bdw}{bandwidth}

\item{kern}{kernel function to use to compute the local-linear weights. Two alternatives: either "gauss" for Gaussian
kernel or "t" for triweight kernel Default: "gauss".}

\item{n1grid}{number of points in the prediction grid for covariate 1. Default at 30.}

\item{n2grid}{number of points in the prediction grid for covariate 2. Default at 30.}
}
\value{
tau matrix ngrid x 3 with median and posterior credible intervals of level 0.90 for the conditional
Kendall's tau.

x1grid grid of n1grid values relative to covariate 1 to each row of the matrix tau

x2grid grid of n2grid values relative to covariate 2 to each row of the matrix tau
}
\description{
This function provides a nonparametric approximation of the posterior distribution of the conditional
Kendall's tau in presence of two covariates.
This function produces estimates using the inconsistent estimator of the conditional copula.
}
\keyword{CondCop}
