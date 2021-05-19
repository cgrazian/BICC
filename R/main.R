#' Nonparametric Bayesian estimate of the conditional Spearman's rho through the empirical likelihood. 
#'
#' This function provides a nonparametric approximation of the posterior distribution of the conditional Spearman's rho. 
#' This function produces estimates using the inconsistent estimator of the conditional copula. 
#'
#' @param uu n x 2 matrix of pseudo-observations. n is the sample size.
#' @param x vector of n values of the one-dimensional covariate
#' @param meth weights to use in the computation of the conditional copula. Two alternatives: "NW" for 
#' Nadaraya-Watson weights, "LL" for local-linear weights
#' @param nsim number of Monte Carlo simulations. Default at 10,000.
#' @param bdw bandwidth
#' @param kern kernel function to use to compute the local-linear weights. Two alternatives: either "gauss" for Gaussian
#' kernel or "t" for triweight kernel Default: "gauss".
#' @param ngrid number of points in the prediction grid for the covariate. Default at 30. 
#' @return rho matrix ngrid x 3 with median and posterior credible intervals of level 0.90 for the conditional 
#' Spearman's rho. 
#' @return x.grid grid of ngrid values relative to each row of the matrix rho 
#' @keywords CondCop
#' @import mnormt
#' @export
ELCOP.rhocov=function(uu,x,meth,nsim=10000,bdw,kern="gauss",ngrid=30)
{
  
  #	install.packages("tmvtnorm")
  n=dim(uu)[1] 
  p=dim(uu)[2]
  ldet.K=-Inf
  
  x.grid=seq( from= (min(x)-sd(x)) , to=(max(x)+sd(x)) , length=ngrid )
  
  K.matern=varcov.Matern(x.grid)
  if(det(K.matern)<0)
  {
    ldet.K=-Inf
  } else {
    ldet.K=log(det(K.matern))
  }
  
  rh.mat=matrix(NA,nrow=nsim,ncol=length(x.grid))
  omega.mat=matrix(NA,nrow=nsim,ncol=length(x.grid))
  
  for(i in 1:nsim)
  {
    
    # Simulate rho from the Gaussian process prior
    w.norm=rmnorm(n = 1, mean = rep(0, ngrid), varcov=K.matern) 
    rh.mat[i,]=(exp(w.norm)-1) / (exp(w.norm)+1)	
    
    # Compute the estimates of rho for the fix x
    for(count in 1:ngrid)
    {
      estim=rho.cond(uu=uu,meth=meth,x.vec=x,xgrid=x.grid[count],
                     bdw=bdw,kern=kern) - rh.mat[i,count]
      omega.mat[i,count]=exp(-EL(estim)$elr)			
    }
    
    print(i)
  }		
  
  rho.cov=matrix(NA,nrow=ngrid,ncol=3)
  for(j in 1:ngrid)
  {
    psam=sample(rh.mat[,j], size=nsim, rep=T, prob=omega.mat[,j])
    summ=c(quantile(psam,.05), quantile(psam,.5),  quantile(psam,.95))
    rho.cov[j,]=sum
  }
  
  return(list(rho=rho.cov,x=x.grid))  
}

#' Nonparametric Bayesian estimate of the conditional Kendall's tau through the empirical likelihood. 
#'
#' This function provides a nonparametric approximation of the posterior distribution of the conditional 
#' Kendall's tau. 
#' This function produces estimates using the inconsistent estimator of the conditional copula. 
#'
#' @param uu n x 2 matrix of pseudo-observations. n is the sample size.
#' @param x vector of n values of the one-dimensional covariate
#' @param meth weights to use in the computation of the conditional copula. Two alternatives: "NW" for 
#' Nadaraya-Watson weights, "LL" for local-linear weights
#' @param nsim number of Monte Carlo simulations. Default at 10,000.
#' @param bdw bandwidth
#' @param kern kernel function to use to compute the local-linear weights. Two alternatives: either "gauss" for Gaussian
#' kernel or "t" for triweight kernel Default: "gauss".
#' @param ngrid number of points in the prediction grid for the covariate. Default at 30. 
#' @return tau matrix ngrid x 3 with median and posterior credible intervals of level 0.90 for the conditional 
#' Kendall's tau.  
#' @return x.grid grid of ngrid values relative to each row of the matrix tau
#' @keywords CondCop
#' @import mnormt
#' @export
ELCOP.taucov=function(uu,x,meth,nsim=10000,bdw,kern="gauss",alpha.gp=0.1,ngrid=30)
{
  
  n=dim(uu)[1] 
  p=dim(uu)[2]
  ldet.K=-Inf
  
  x.grid=seq( from= (min(x)-sd(x)) , to=(max(x)+sd(x)) , length=ngrid )
  
  K.matern=varcov.Matern(x.grid)
  if(det(K.matern)<0)
  {
    ldet.K=-Inf
  } else {
    ldet.K=log(det(K.matern))
  }
  
  tau.mat=matrix(NA,nrow=nsim,ncol=length(x.grid))
  omega.mat=matrix(NA,nrow=nsim,ncol=length(x.grid))
  
  for(i in 1:nsim)
  {
    
    # Simulate rho from the Gaussian process prior
    w.norm=rmnorm(n = 1, mean = rep(0, ngrid), varcov=K.matern) 
    tau.mat[i,]=(exp(w.norm)-1) / (exp(w.norm)+1)	
    
    # Compute the estimates of rho for the fix x
    for(count in 1:ngrid)
    {
      estim=tau.cond(uu=uu,meth=meth,x.vec=x,xgrid=x.grid[count],
                     bdw=bdw,kern=kern) - tau.mat[i,count]
      omega.mat[i,count]=exp(-EL(estim)$elr)			
    }
    
    print(i)
  }		
  
  tau.cov=matrix(NA,nrow=ngrid,ncol=3)
  for(j in 1:ngrid)
  {
    psam=sample(tau.mat[,j], size=nsim, rep=T, prob=omega.mat[,j])
    sintesi=c(quantile(psam,.05), quantile(psam,.5),  quantile(psam,.95))
    tau.cov[j,]=sintesi
  }
  
  return(list(tau=tau.cov,x=x.grid))  
}


#' Nonparametric Bayesian estimate of the conditional Spearman's rho through the empirical likelihood, with 
#' a linearised model. 
#'
#' This function provides a nonparametric approximation of the posterior distribution of the conditional Spearman's rho. 
#' This function produces estimates computing the likelihood for the linearised model of the conditional Spearman's rho, 
#' obtained through a Taylor expansion. 
#'
#' @param w input data, i.e. Fisher transform of the estimated Spearman's rho. 
#' @param x vector of n values of the one-dimensional covariate
#' @param p order of the approximation in the Taylor expansion.
#' @param nsim number of Monte Carlo simulations. Default at 10,000.
#' @param mu prior mean of the coefficients. If NULL, the mean will be a vector of zeros of length p. Default at NULL
#' @param Sigma prior covariance matrix of the coefficients of the linearised model. If NULL, it is fixed at a diagonal
#' matrix with elements on the diagonal equal to 10. Default at NULL.
#' @return beta: a matrix of dimension nsim x p with a sample of coefficients parameters
#' @return omega: a vector of nsim weights representing the empirical likelihood for the coefficients in the matrix beta. 
#' @keywords CondCop
#' @import mvtnorm
#' @export
ELCOP.cond=function(w,x,p,nsim=10000,mu=NULL,Sigma=NULL)
{
  n <- length(w)
  x.des <- matrix(NA, ncol=(p+1), nrow=n)
  x.des[,1] <- rep(1,n)
  
  for(d in 2:(p+1))
  {
    x.des[,d] <- x^(d-1)
  }
  
  if(is.null(mu)==T){
    mu <- rep(0,(p+1))
  }
  
  if(is.null(Sigma)==T){
    Sigma <- 10 * diag((p+1))
  }
  
  beta_mat=matrix(NA,nrow=nsim,ncol=(p+1))
  omega.mat <- matrix(NA,nsim,1)
  
  for(j in 1:nsim)
  {
    
    beta_prop <- rmvnorm(1,mu, Sigma)
    beta_mat[j,] <- beta_prop    
    beta_prop <- t(beta_prop)
    
    estim <- matrix(0,nrow=n,ncol=(p+1))
    
    for(i in 1:n)
    {
      x.i <- matrix(x.des[i,],ncol=1,nrow=(p+1))
      estim[i,] <- x.i %*% (w[i] - t(x.i) %*% beta_prop)
    }
    cond <- apply(estim,2,mean)
    
    omega.mat[j,] <- exp(-EL(cond)$elr)
    
    if(round(j/10000)==j/10000){
      print(j)
    }
  }		
  
  return(list(beta=beta_mat, omega=omega.mat))  
}

#' Nonparametric Bayesian estimate of the conditional Kendall's tau through the empirical likelihood. 
#'
#' This function provides a nonparametric approximation of the posterior distribution of the conditional 
#' Kendall's tau in presence of two covariates.
#' This function produces estimates using the inconsistent estimator of the conditional copula. 
#'
#' @param uu n x 2 matrix of pseudo-observations. n is the sample size.
#' @param x1 vector of n values of covariate 1
#' @param x2 vector of n values of covariate 2
#' @param nsim number of Monte Carlo simulations. Default at 10,000.
#' @param bdw bandwidth
#' @param kern kernel function to use to compute the local-linear weights. Two alternatives: either "gauss" for Gaussian
#' kernel or "t" for triweight kernel Default: "gauss".
#' @param n1grid number of points in the prediction grid for covariate 1. Default at 30. 
#' @param n2grid number of points in the prediction grid for covariate 2. Default at 30. 
#' @return tau matrix ngrid x 3 with median and posterior credible intervals of level 0.90 for the conditional 
#' Kendall's tau.
#' @return x1grid grid of n1grid values relative to covariate 1 to each row of the matrix tau 
#' @return x2grid grid of n2grid values relative to covariate 2 to each row of the matrix tau 
#' @keywords CondCop
#' @import mnormt
#' @export
ELCOP.taucov_multi=function(uu,x1,x2,nsim=10000,bdw,kern="gauss",
                            n1grid=30,n2grid=30)
{
  
  n=dim(uu)[1] 
  p=dim(uu)[2]
  
  x1.grid=seq( from= (min(x1)-sd(x1)) , to=(max(x1)+sd(x1)) , length=n1grid )
  x2.grid=seq( from= (min(x2)-sd(x2)) , to=(max(x2)+sd(x2)) , length=n2grid )
  
  tau.list <- list()
  omega.list <- list()
  
  for(i in 1:nsim)
  {
    
    # Simulate rho from the Gaussian process prior
    w.norm=matrix(rnorm(n1grid*n2grid,0,1),n1grid,n2grid) 
    tau.mat=(exp(w.norm)-1) / (exp(w.norm)+1)	
    
    # Compute the estimates of rho for the fix x
    tau_estim <- tau.cond_multi(uu=uu,cbind(x1,x2),x1.grid,x2.grid,bdw,kern)
    omega_mat <- matrix(NA,ncol=ncol(tau_estim),nrow=nrow(tau_estim))
    for(d1 in 1:nrow(tau_estim))
    {
      for(d2 in 1:ncol(tau_estim))
      {
        estim=tau_estim[d1,d2] - tau.mat[d1,d2]
        omega_mat[d1,d2]=exp(-EL(estim)$elr)			
      }
    }
    tau.list[[i]] <- tau.mat
    omega.list[[i]] <- omega_mat
  }		
  
  return(list(tau=tau.list,omega=omega.list,x1grid=x1.grid,x2grid=x2.grid))  
}
