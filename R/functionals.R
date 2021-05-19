#' Conditional Spearman's rho.
#'
#' This function provides a nonparametric point estimate for the conditional Spearman'rho
#'
#' @param uu n x 2 matrix of pseudo-observations. n is the sample size.
#' @param meth weights to use in the computation of the conditional copula. Two alternatives: "NW" for
#' Nadaraya-Watson weights, "LL" for local-linear weights
#' @param x.vec vector of n values of the one-dimensional covariate
#' @param xgrid vector representing a grid of values with respect to which to compute the weights
#' @param bdw bandwidth
#' @param kern kernel function to use to compute the local-linear weights. Two alternatives: either "gauss" for Gaussian
#' kernel or "t" for triweight kernel Default: "gauss".
#' @return nonparametric point estimate of the Spearman's rho.
#' @keywords CondCop
#' @export
rho.cond=function(uu,meth="NW",x.vec,xgrid,bdw,kern="gauss")
{
  uu1=uu[,1]
  uu2=uu[,2]
  estim=c()

  for(count in 1:length(xgrid))
  {
    # Weights
    if(meth=="NW")
    {
      ww=NW.weights(x.vec,xgrid[count],bdw,kern=kern)
    } else {
      ww=LL.weights(x.vec,xgrid[count],bdw,kern=kern)
    }

    # Compute uu.hat
    uu1.hat=c()
    uu2.hat=c()

    for(i in 1:nrow(uu))
    {
      uu1.hat[i]=sum(ww*(uu1<=uu1[i]))
      uu2.hat[i]=sum(ww*(uu2<=uu2[i]))
    }

    estim[count]=12*sum(ww*(1-uu1.hat)*(1-uu2.hat))-3
  }

  # Return
  return(estim)
}

#' Conditional Kendall's tau.
#'
#' This function provides a nonparametric point estimate for the conditional Kendall's tau
#'
#' @param uu n x 2 matrix of pseudo-observations. n is the sample size.
#' @param meth weights to use in the computation of the conditional copula. Two alternatives: "NW" for
#' Nadaraya-Watson weights, "LL" for local-linear weights
#' @param x.vec vector of n values of the one-dimensional covariate
#' @param xgrid vector representing a grid of values with respect to which to compute the weights
#' @param bdw bandwidth
#' @param kern kernel function to use to compute the local-linear weights. Two alternatives: either "gauss" for Gaussian
#' kernel or "t" for triweight kernel Default: "gauss".
#' @return nonparametric point estimate of the Kendall's tau.
#' @keywords CondCop
#' @export
tau.cond=function(uu,meth="NW",x.vec,xgrid,bdw,kern="gauss")
{
  uu1=uu[,1]
  uu2=uu[,2]
  estim=c()

  for(count in 1:length(xgrid))
  {
    # Weights
    if(meth=="NW")
    {
      ww=NW.weights(x.vec,xgrid[count],bdw,kern=kern)
    } else {
      ww=LL.weights(x.vec,xgrid[count],bdw,kern=kern)
    }

    # Compute uu.hat
    uu1.hat=c()
    uu2.hat=c()

    for(i in 1:nrow(uu))
    {
      uu1.hat[i]=sum(ww*(uu1<=uu1[i]))
      uu2.hat[i]=sum(ww*(uu2<=uu2[i]))
    }

    estim[count]=12*sum(ww*(1-uu1.hat)*(1-uu2.hat))-3
  }

  # Return
  return(estim)
}

#' Conditional Kendall's tau for multivariate covariates.
#'
#' This function provides a nonparametric point estimate for the conditional Kendall's tau
#'
#' @param uu n x 2 matrix of pseudo-observations. n is the sample size.
#' @param meth weights to use in the computation of the conditional copula. Two alternatives: "NW" for
#' Nadaraya-Watson weights, "LL" for local-linear weights
#' @param x.mat nx2 matrix of values of a 2-dimensional covariate
#' @param x1.grid vector representing a grid of values for covariate 1
#' @param x2.grid vector representing a grid of values for covariate 2
#' @param bdw bandwidth
#' @param kern kernel function to use to compute the NW weights. Two alternatives: either "gauss" for Gaussian
#' kernel or "t" for triweight kernel Default: "gauss".
#' @return matrix with nonparametric point estimate of the Kendall's tau for each combination of the values
#' in x1.grid and x2.grid
#' @keywords CondCop
#' @export
tau.cond_multi=function(uu,x.mat,x.grid,bdw,kern="gauss")
{
  uu1=uu[,1]
  uu2=uu[,2]
  tau_mat=matrix(NA,nrow=length(x1.grid),ncol=length(x2.grid))


  for(d1 in 1:length(x1.grid))
  {
    for(d2 in 1:length(x2.grid))
    {
      # Weights
      ww=NW.weights_multi(x.mat,c(x1.grid[d1],x2.grid[d2]),bdw,kern=kern)
      ww <- ifelse(is.na(ww)==T,0,ww)

      # Compute Kendall's tau
      sum_cop <- 0
      for(i in 1:nrow(uu)){
        for(j in 1:nrow(uu)){
          sum_cop <- sum_cop + ww[i]*ww[j]* (uu1[i]<uu[j] & uu2[i]<uu2[j])
        }
      }
      tau_mat[d1,d2]= -1 + 4/(1-sum(ww^2)) * sum_cop
    }
  }

  # Return
  return(tau_mat)
}

