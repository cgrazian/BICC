#' Quadratic function for the empirical likelihood
#'
#' This function implements the quadratic function in Owen (2001) (Section 3.14).
#' This is the function to optimise to obtain the empirical likelihood. 
#'
#' @param x observations
#' @param thre threshold 
#' @param gradient gradient
#' @return A vector representing the pseudo-logaritm of the values to minimise in the empirical likelihood
#' @keywords CondCop
#' @export

Owen <- function(x, thre, gradient)
{
  grad <- t(gradient)
  par <- t(as.matrix(x))
  eps <- 1/thre
  
  z <- 1+par%*%grad
  ans <- z
  lo <- ((z)<eps)
  
  ans[ lo  ] <- log(eps) - 1.5 + 2*(z[lo])/eps - 0.5*(z[lo]/eps)^2
  ans[ !lo ] <- log( z[!lo] )
  -sum(ans)
}

#' Empirical likelihood
#'
#' This function computes the empirical likelihood by minimising a quadratic function. 
#'
#' @param SC individual contributions of the estimating function. 
#' It can be either an n x 1 vector (scalar parameter) or an n x p matrix (vector valued parameter)
#' @return w0: weights associated to each unit
#' @return conv: convergence of the algorithm (available only for vector valued parameter)
#' @return elr: empirical loglikelihood ratio
#' @return el: empirical loglikelihood
#' @keywords CondCop
#' @export
EL <- function(SC)
{
  n <- NROW(SC)
  p <- NCOL(SC)
  
  ##find Lagrange multiplier
  if(p!=1)
  {
    OBJ <- optim(par=rep(0,p), fn=Owen, thre=n, gradiente=SC, control=list(maxit=1e4))
    molt_lagr_0 <- OBJ$pa
    conv <- OBJ$conv
  }
  else
  {
    molt_lagr_0 <- optim(par=rep(0,p), fn=Owen, thre=n, gradiente=SC, method="Brent", lower=-1e1, upper=1e1)$pa
  }
  
  ##weights
  w_emp_0 <- as.numeric( 1/(n*(1+molt_lagr_0%*%t(SC))) )
  
  if(p!=1)
  {
    list(w0=w_emp_0, conv=conv, elr=-2*sum(log(w_emp_0))-2*n*log(n), el=log(w_emp_0))
  }
  else
  {
    list(w0=w_emp_0, elr=-2*sum(log(w_emp_0))-2*n*log(n), el=log(w_emp_0) )
  }
}

#' Triweight kernel function
#'
#' This function computes the triweight kernel function of an input vector
#'
#' @param x scalar or vector on which to compute the triweight function
#' @return value of the triweight function
#' @keywords CondCop
#' @export
triweight=function(x)
{
  35/32*(1-x^2)^3 *(abs(x)<1)
}

#' Gaussian kernel function
#'
#' This function computes the Gaussian kernel function of an input vector
#'
#' @param x scalar or vector on which to compute the Gaussian function
#' @return value of the Gaussian function
#' @keywords CondCop
#' @export
gaussweight=function(x)
{
  exp(-x^2/2)/sqrt(2*pi)
}

#' Nadaraya-Watson weights
#'
#' This function computes the Nadaraya-Watson function weights 
#'
#' @param x.vec vector on which to compute the Nadaraya-Watson weights
#' @param x scalar: this is the point with respect to which to compute the Nadaraya-Watson weights for each of the points of x.vec
#' @param bdw bandwidth
#' @param kern kernel function to use to compute the Nadaraya-Watson weights. Two alternatives: either "gauss" for Gaussian
#' kernel or "t" for triweight kernel Default: "gauss".
#' @return vector of weights
#' @keywords CondCop
#' @export
NW.weights=function(x.vec,x,bdw,kern="gauss")
{
  val=(x.vec-x)/bdw
  if(kern=="t"){
    const=sum(triweight(val))
    temp=triweight(val)/const
  } else {
    const=sum(gaussweight(val))
    temp=gaussweight(val)/const
  }
  return(temp)
}

#' Multivariate Nadaraya-Watson weights
#'
#' This function computes the multivariate Nadaraya-Watson function weights 
#'
#' @param x.mat matrix of value on which to compute the Nadaraya-Watson weights
#' @param x.vec vector: this is the point with respect to which to compute the Nadaraya-Watson weights for each of the points of x.vec
#' @param bdw bandwidth
#' @param kern kernel function to use to compute the Nadaraya-Watson weights. Two alternatives: either "gauss" for Gaussian
#' kernel or "t" for triweight kernel Default: "gauss".
#' @return vector of weights
#' @keywords CondCop
#' @export
NW.weights_multi=function(x.mat,x.vec,bdw,kern="gauss")
{
  d <- length(x.vec)
  n <- nrow(x.mat)
  NW_mat <- matrix(NA,nrow=n,ncol=d) 
  for(di in 1:d)
  {
    val <- (x.mat[,di] - x.vec[di]) / bdw   
    if(kern=="t"){
      temp=triweight(val)
    } else {
      temp=gaussweight(val)
    }
    NW_mat[,di] <- temp
  }
  lNW_mat <- log(NW_mat)
  sum_lNW_mat <- apply(lNW_mat,1,sum)
  #    NW_vec <- apply(NW_mat,1,prod)
  NW_vec <- exp(sum_lNW_mat)/sum(exp(sum_lNW_mat))
  return(NW_vec)
}

#' Local-linear weights
#'
#' This function computes the local-linear function weights 
#'
#' @param x.vec vector on which to compute the local-linear weights
#' @param x scalar: this is the point with respect to which to compute the local-linear weights for each of the 
#' points of x.vec
#' @param bdw bandwidth
#' @param kern kernel function to use to compute the local-linear weights. Two alternatives: either "gauss" for Gaussian
#' kernel or "t" for triweight kernel Default: "gauss".
#' @return vector of weights
#' @keywords CondCop
#' @export
LL.weights=function(x.vec,x,bdw,kern="gauss")
{
  n=length(x.vec)
  val=(x.vec-x)/bdw
  
  if(kern=="t"){
    K=triweight(val)
  } else {
    K=gaussweight(val)
  }
  
  S0=sum((val^0)*K)/(n*bdw)
  S1=sum((val^1)*K)/(n*bdw)
  S2=sum((val^2)*K)/(n*bdw)
  temp=(K*(S2-val*S1))/(n*bdw*(S0*S2-S1^2))
  
  return(temp)
}

#' Matern covariance function
#'
#' This function computes the Matern covariance function for a vector X.
#'
#' @param X vector to compute the covariance matrix 
#' @param nu scalar; smoothness parameter of the Matern covariance matrix. Default at 5/2. 
#' @return a matrix of the dimension n x n where n is the length of X. 
#' @keywords CondCop
#' @import fields 
#' @export
varcov.Matern=function(X,nu=5/2){
  K=matrix(nrow=length(X),ncol=length(X))
  for(i in 1:length(X)){
    for(j in 1:length(X)){
      d=sqrt((X[i]-X[j])^2)
      K[i,j]=Matern(d,nu=nu)
    }
  }
  return(K)
}

#' Fisher transform
#'
#' This function computes the Fisher transform for a vector y.
#'
#' @param y vector to compute the Fisher transform
#' @return a vector of the same length of y with the values of the Fisher transform
#' @keywords CondCop
#' @import fields 
#' @export
fish_transf <- function(y){
  0.5 * log( (1+y) / (1-y) )
}
