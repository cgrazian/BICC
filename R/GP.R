#' Gaussian process modelling for the Fisher transform of the condtional functionals.
#'
#' This function makes use of the tgp package developed by Robert B. Gramacy (Virginia Polytechnic and State University)
#' and Dr. Matt A. Taddy (Amazon)
#' https://cran.r-project.org/web/packages/tgp/tgp.pdf
#'
#'X=z, XX=x.grid, Z=W_rep
#' @param z values of the covariate
#' @param nz number of values in the grid of values in the domain of the covariate. Default at 1000. 
#' @param y response variable. Typically this is the Fisher transform of the estimated conditional functional. 
#' @param model Bayesian regression models to be implemented. "blm" stands for linear model, 
#' "bgpllm" stands for GP with jumps to the limiting linear model (LLM), 
#' "bgp" stands for Gaussian process. Default at "bgp"
#' @return a list including 
#' estim: point estimate for the functional of interest phi
#' q1: quantile of level 0.025 for the functional of interest phi
#' q2: quantile of level 0.975 for the functional of interest phi
#' xgrid: a vector of points representing a grid of values on which to approximate the predictive distribution
#' @return x.grid grid of ngrid values relative to each row of the matrix rho 
#' @keywords CondCop
#' @import tgp
#' @export
cond_GP <- function(z,nz,y, model="bgp"){
  x.grid=seq( from= (min(z)-sd(z)) , to=(max(z)+sd(z)) , length=1000 )
  if(model=="blm"){
    obj <- blm(X=z, XX=x.grid, Z=W_rep)    # too linear
  } else {
    if(model=="bgpllm"){
      obj <- bgpllm(X=z, XX=x.grid, Z=W_rep)    
    } else {
      obj <- bgp(X=z, Z=W_rep, XX=x.grid, verb=0)
    }
  }
  
  phi_GP <- (exp(2*obj$ZZ.med)-1) / (exp(2*obj$ZZ.med)+1)
  phi_GP.q1 <- (exp(2*obj$ZZ.q1)-1) / (exp(2*obj$ZZ.q1)+1)
  phi_GP.q2 <- (exp(2*obj$ZZ.q2)-1) / (exp(2*obj$ZZ.q2)+1) 
  
  return(list(estim=phi_GP,q1=phi_GP.q1,q2=phi_GP.q2,xgrid=x.grid))
}
