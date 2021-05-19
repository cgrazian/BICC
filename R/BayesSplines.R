#' Bayesian fitting of shape-restricted regression spline.
#'
#' This function fit a Bayesian shape-restricted regression spline, using a user-chosen algorithm type.
#' Code available on Prof. Mary C. Meyer's website (Colorado State University) https://www.stat.colostate.edu/~meyer/bayescode.htm
#' @param type type of the algorithm (and the restriction on the shape)
#' 1 = monotone increasing (quadratic)
#  2 = monotone decreasing (quadratic)
#  3 = convex (cubic)
#  4 = concave (cubic)
#  5 = increasing convex (cubic)
#  6 = decreasing convex (cubic)
#  7 = increasing concave (cubic)
#  8 = decreasing concave (cubic)
#' @param x vector of n values of the one-dimensional covariate
#' @param y vector of n values of the response variable
#' @param k number of interior knots
#' @param nloop number of loops in the algorithm
#' @return  a list of elements
#' fit: point fitted values
#' lower: lower bound of the credible interval
#' upper: upper bound of the credible interval
#' coefs: coefficients of the spline regression
#' @keywords CondCop
#' @export
brspl=function(type,x,y,k,nloop){
  sdy=sqrt(var(y))
  n=length(x)
  xs=sort(x)
  knots=round((0:(k+1))*n/(k+1))
  knots[1]=1
  t=xs[knots]
  if(type==1){delta=monincr(x,t)}
  if(type==2){delta=mondecr(x,t)}
  if(type==3){delta=convex(x,t)}
  if(type==4){delta=concave(x,t)}
  if(type==5){delta=incconvex(x,t)}
  if(type==7){delta=incconcave(x,t)}
  if(type==6){delta=decconvex(x,t)}
  if(type==8){delta=decconcave(x,t)}
  incr=0
  decr=0
  if(type==1|type==5|type==7){incr=1}
  if(type==2|type==6|type==8){decr=1}
  m=length(delta)/n
  s1=1:m
  for(i in 1:m){s1[i]=sum(delta[i,]^2)}
  v=1:n*0+1
  vx=x-mean(x)
  sumvx=sum(vx^2)
  xmat=cbind(v,vx)
  line=xmat%*%solve(t(xmat)%*%xmat)%*%t(xmat)%*%y
  xmatb=cbind(v,t(delta))
  if(incr==0&decr==0){xmatb=cbind(xmatb,vx)}
  range=sdy
  # initialize
  #nloop=5000
  beta=matrix(1:(m*nloop)*0,nrow=m)
  tau=1:nloop*0
  alpha=1:nloop*0
  if(decr==0&incr==0){alphax=1:nloop*0;alphax[1]=0}
  beta[,1]=1:m*0+range/m
  tau[1]=1
  alpha[1]=mean(y)
  r=1:n*0
  # Set prior parameters
  #ssq=sseb/(n-kb)
  d1=1/10
  d2=1/10
  c1=1/m^2
  a=c1-1
  c2=1/m/range
  mvar=1000
  #loop!
  for(iloop in 2:nloop){
    #first get a new alpha:
    r=y-t(delta)%*%beta[,iloop-1]
    if(decr==0&incr==0){r=r-alphax[iloop-1]*vx}
    amean=tau[iloop-1]*sum(r)/(1/mvar+n*tau[iloop-1])
    astd=(1/mvar+n*tau[iloop-1])^(-1/2)
    alpha[iloop]=rnorm(1,amean,astd)
    #for convex and concave, need new alphax
    if(decr==0&incr==0){
      r=y-t(delta)%*%beta[,iloop-1]-alpha[iloop]
      s2=sum(r*vx)
      amean=tau[iloop-1]*s2/(1/mvar+tau[iloop-1]*sumvx)
      astd=(1/mvar+tau[iloop-1]*sumvx)^(-1/2)
      alphax[iloop]=rnorm(1,amean,astd)
    }
    #next find a new tau:
    r=y-t(delta)%*%beta[,iloop-1]-alpha[iloop]
    if(incr==0&decr==0){r=r-alphax[iloop]*vx}
    sse=sum(r^2)
    nd1=d1+n/2
    nd2=d2+sse/2
    tau[iloop]=rgamma(1,nd1,nd2)
    #  loop through the betas: (random order)
    list=order(runif(m))
    for(j0 in 1:m){
      j=list[j0]
      index=list!=j
      buse=beta[list[index],iloop-1]
      buse2=beta[list[index],iloop]
      if(j0>1){
        for(l in 1:(j0-1)){buse[l]=buse2[l]}
      }
      r=y-alpha[iloop]-t(delta[list[index],])%*%buse
      if(decr==0&incr==0){r=r-alphax[iloop]*vx}
      s2=sum(r*delta[j,])
      b=tau[iloop]*s1[j]/2
      c=s2/s1[j]-c2/tau[iloop]/s1[j]
      beta[j,iloop]=postgen(a,b,c)
    }
  }
  beta0=1:m
  warm=round(nloop/10)
  use=nloop-warm
  for(j in 1:m){beta0[j]=mean(beta[j,warm:nloop])}
  alpha0=mean(alpha)
  bfit=(alpha0+t(delta)%*%beta0)
  if(incr==0&decr==0){bfit=bfit+mean(alphax)*vx}
  # make credible intervals
  fits=matrix(nrow=n,ncol=use)
  for(j in 1:use){
    i=warm+j
    fits[,j]=t(delta)%*%beta[,i]+alpha[i]
    if(incr==0&decr==0){fits[,j]=fits[,j]+alphax[i]*vx}
  }
  f1=1:n;f2=1:n;fm=1:n
  blower=round(use*.025)
  bupper=round(use*.975)
  bmid=round(use/2)
  for(i in 1:n){
    fs=sort(fits[i,])
    f1[i]=fs[blower]
    f2[i]=fs[bupper]
    fm[i]=fs[bmid]
  }
  ans=new.env()
  ans$fit=bfit
  ans$lower=f1
  ans$upper=f2
  ans$coefs=beta0
  ans
}


#' Posterior generation for the beta coefficients in shape-restricted splines.
#'
#' @param a parameter a in the regression spline model
#' @param b parameter b in the regression spline model
#' @param c parameter c in the regression spline model
#' @return transformed value for beta
#' @keywords CondCop
#' @export
postgen=function(a,b,c){
  x=c
  if(c<1e-8){x=1e-8}
  a1=a+1
  for(i in 1:5){
    u=runif(1)*exp(-b*(x-c)^2)
    bl=max(0,c-sqrt(-log(u)/b))
    bu=c+sqrt(-log(u)/b)
    e=runif(1)
    x=(e*bu^a1+(1-e)*bl^a1)^(1/a1)
  }
  x
}

#' Monotone increasing function
#'
#' @param x value of the covariate
#' @param t value of t
#' @return value
#' @keywords CondCop
#' @export
monincr=function(x,t){
  n=length(x)
  k=length(t)-2
  m=k+2
  sigma=matrix(1:m*n,nrow=m,ncol=n)
  for(j in 1:(k-1)){
    i1=x<=t[j]
    sigma[j,i1] = 0
    i2=x>t[j]&x<=t[j+1]
    sigma[j,i2] = (x[i2]-t[j])^2 / (t[j+2]-t[j]) / (t[j+1]-t[j])
    i3=x>t[j+1]&x<=t[j+2]
    sigma[j,i3] = 1-(x[i3]-t[j+2])^2/(t[j+2]-t[j+1])/(t[j+2]-t[j])
    i4=x>t[j+2]
    sigma[j,i4]=1
  }

  i1=x<=t[k]
  sigma[k,i1] = 0
  i2=x>t[k]&x<=t[k+1]
  sigma[k,i2] = (x[i2]-t[k])^2 / (t[k+2]-t[k]) / (t[k+1]-t[k])
  i3=x>t[k+1]&x<=t[k+2]
  sigma[k,i3] = 1- (x[i3]-t[k+2])^2/(t[k+2]-t[k+1])/(t[k+2]-t[k])
  i4=x>t[k+2]
  sigma[k,i4]=1

  i1=x<=t[2]
  sigma[k+1,i1]=1-(t[2]-x[i1])^2/(t[2]-t[1])^2
  i2=x>t[2]
  sigma[k+1,i2]=1

  i1=x<=t[k+1]
  sigma[k+2,i1]=0
  i2=x>t[k+1]&x<=t[k+2]
  sigma[k+2,i2]=(x[i2]-t[k+1])^2/(t[k+2]-t[k+1])^2
  i3=x>t[k+2]
  sigma[k+2,i3]=1
  for(i in 1:m){sigma[i,]=sigma[i,]-mean(sigma[i,])}
  sigma
}

#' Monotone decreasing function
#'
#' @param x value of the covariate
#' @param t value of t
#' @return value
#' @keywords CondCop
#' @export
mondecr=function(x,t){
  n=length(x)
  k=length(t)-2
  m=k+2
  sigma=matrix(1:m*n,nrow=m,ncol=n)
  for(j in 1:(k-1)){
    i1=x<=t[j]
    sigma[j,i1] = 0
    i2=x>t[j]&x<=t[j+1]
    sigma[j,i2] = (x[i2]-t[j])^2 / (t[j+2]-t[j]) / (t[j+1]-t[j])
    i3=x>t[j+1]&x<=t[j+2]
    sigma[j,i3] = 1-(x[i3]-t[j+2])^2/(t[j+2]-t[j+1])/(t[j+2]-t[j])
    i4=x>t[j+2]
    sigma[j,i4]=1
  }

  i1=x<=t[k]
  sigma[k,i1] = 0
  i2=x>t[k]&x<=t[k+1]
  sigma[k,i2] = (x[i2]-t[k])^2 / (t[k+2]-t[k]) / (t[k+1]-t[k])
  i3=x>t[k+1]&x<=t[k+2]
  sigma[k,i3] = 1- (x[i3]-t[k+2])^2/(t[k+2]-t[k+1])/(t[k+2]-t[k])
  i4=x>t[k+2]
  sigma[k,i4]=1

  i1=x<=t[2]
  sigma[k+1,i1]=1-(t[2]-x[i1])^2/(t[2]-t[1])^2
  i2=x>t[2]
  sigma[k+1,i2]=1

  i1=x<=t[k+1]
  sigma[k+2,i1]=0
  i2=x>t[k+1]&x<=t[k+2]
  sigma[k+2,i2]=(x[i2]-t[k+1])^2/(t[k+2]-t[k+1])^2
  i3=x>t[k+2]
  sigma[k+2,i3]=1
  for(i in 1:m){sigma[i,]=sigma[i,]-mean(sigma[i,])}
  sigma=-sigma
  sigma
}


#' Convex function
#'
#' @param x value of the covariate
#' @param t value of t
#' @return value
#' @keywords CondCop
#' @export
convex=function(x,t){
  n=length(x)
  k=length(t)-2
  m=k+2
  sigma=matrix(1:m*n,nrow=m,ncol=n)
  for(j in 1:(k-1)){
    i1=x<=t[j]
    sigma[j,i1] = 0
    i2=x>t[j]&x<=t[j+1]
    sigma[j,i2] = (x[i2]-t[j])^3 / (t[j+2]-t[j]) / (t[j+1]-t[j])/3
    i3=x>t[j+1]&x<=t[j+2]
    sigma[j,i3] = x[i3]-t[j+1]-(x[i3]-t[j+2])^3/(t[j+2]-t[j])/(t[j+2]-t[j+1])/3+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
    i4=x>t[j+2]
    sigma[j,i4]=(x[i4]-t[j+1])+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
  }
  i1=x<=t[k]
  sigma[k,i1] = 0
  i2=x>t[k]&x<=t[k+1]
  sigma[k,i2] = (x[i2]-t[k])^3 / (t[k+2]-t[k]) / (t[k+1]-t[k])/3
  i3=x>t[k+1]
  sigma[k,i3] = x[i3]-t[k+1]-(x[i3]-t[k+2])^3/(t[k+2]-t[k])/(t[k+2]-t[k+1])/3+(t[k+1]-t[k])^2/3/(t[k+2]-t[k])-(t[k+2]-t[k+1])^2/3/(t[k+2]-t[k])
  i1=x<=t[2]
  sigma[k+1,i1]=x[i1]-t[1]+(t[2]-x[i1])^3/(t[2]-t[1])^2/3
  i2=x>t[2]
  sigma[k+1,i2]=x[i2]-t[1]
  i1=x<=t[k+1]
  sigma[k+2,i1]=0
  i2=x>t[k+1]
  sigma[k+2,i2]=(x[i2]-t[k+1])^3/(t[k+2]-t[k+1])^2/3
  v1=1:n*0+1
  v2=x
  xmat=cbind(v1,v2)
  pr=xmat%*%solve(t(xmat)%*%xmat)%*%t(xmat)
  for(j in 1:m){
    sigma[j,]=sigma[j,]-pr%*%sigma[j,]
    sigma[j,]=sigma[j,]/(max(sigma[j,])-min(sigma[j,]))
  }
  sigma
}

#' Concave function
#'
#' @param x value of the covariate
#' @param t value of t
#' @return value
#' @keywords CondCop
#' @export
concave=function(x,t){
  n=length(x)
  k=length(t)-2
  m=k+2
  sigma=matrix(1:m*n,nrow=m,ncol=n)
  for(j in 1:(k-1)){
    i1=x<=t[j]
    sigma[j,i1] = 0
    i2=x>t[j]&x<=t[j+1]
    sigma[j,i2] = (x[i2]-t[j])^3 / (t[j+2]-t[j]) / (t[j+1]-t[j])/3
    i3=x>t[j+1]&x<=t[j+2]
    sigma[j,i3] = x[i3]-t[j+1]-(x[i3]-t[j+2])^3/(t[j+2]-t[j])/(t[j+2]-t[j+1])/3+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
    i4=x>t[j+2]
    sigma[j,i4]=(x[i4]-t[j+1])+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
  }
  i1=x<=t[k]
  sigma[k,i1] = 0
  i2=x>t[k]&x<=t[k+1]
  sigma[k,i2] = (x[i2]-t[k])^3 / (t[k+2]-t[k]) / (t[k+1]-t[k])/3
  i3=x>t[k+1]
  sigma[k,i3] = x[i3]-t[k+1]-(x[i3]-t[k+2])^3/(t[k+2]-t[k])/(t[k+2]-t[k+1])/3+(t[k+1]-t[k])^2/3/(t[k+2]-t[k])-(t[k+2]-t[k+1])^2/3/(t[k+2]-t[k])
  i1=x<=t[2]
  sigma[k+1,i1]=x[i1]-t[1]+(t[2]-x[i1])^3/(t[2]-t[1])^2/3
  i2=x>t[2]
  sigma[k+1,i2]=x[i2]-t[1]
  i1=x<=t[k+1]
  sigma[k+2,i1]=0
  i2=x>t[k+1]
  sigma[k+2,i2]=(x[i2]-t[k+1])^3/(t[k+2]-t[k+1])^2/3
  v1=1:n*0+1
  v2=x
  xmat=cbind(v1,v2)
  pr=xmat%*%solve(t(xmat)%*%xmat)%*%t(xmat)
  for(j in 1:m){
    sigma[j,]=sigma[j,]-pr%*%sigma[j,]
    sigma[j,]=sigma[j,]/(max(sigma[j,])-min(sigma[j,]))
  }
  sigma=-sigma
  sigma
}

#' Increasing convex function
#'
#' @param x value of the covariate
#' @param t value of t
#' @return value
#' @keywords CondCop
#' @export
incconvex=function(x,t){
  n=length(x)
  k=length(t)-2
  m=k+3
  sigma=matrix(1:m*n,nrow=m,ncol=n)
  for(j in 1:(k-1)){
    i1=x<=t[j]
    sigma[j,i1] = 0
    i2=x>t[j]&x<=t[j+1]
    sigma[j,i2] = (x[i2]-t[j])^3 / (t[j+2]-t[j]) / (t[j+1]-t[j])/3
    i3=x>t[j+1]&x<=t[j+2]
    sigma[j,i3] = x[i3]-t[j+1]-(x[i3]-t[j+2])^3/(t[j+2]-t[j])/(t[j+2]-t[j+1])/3+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
    i4=x>t[j+2]
    sigma[j,i4]=(x[i4]-t[j+1])+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
  }
  i1=x<=t[k]
  sigma[k,i1] = 0
  i2=x>t[k]&x<=t[k+1]
  sigma[k,i2] = (x[i2]-t[k])^3 / (t[k+2]-t[k]) / (t[k+1]-t[k])/3
  i3=x>t[k+1]
  sigma[k,i3] = x[i3]-t[k+1]-(x[i3]-t[k+2])^3/(t[k+2]-t[k])/(t[k+2]-t[k+1])/3+(t[k+1]-t[k])^2/3/(t[k+2]-t[k])-(t[k+2]-t[k+1])^2/3/(t[k+2]-t[k])
  i1=x<=t[2]
  sigma[k+1,i1]=x[i1]-t[1]+(t[2]-x[i1])^3/(t[2]-t[1])^2/3
  i2=x>t[2]
  sigma[k+1,i2]=x[i2]-t[1]
  i1=x<=t[k+1]
  sigma[k+2,i1]=0
  i2=x>t[k+1]
  sigma[k+2,i2]=(x[i2]-t[k+1])^3/(t[k+2]-t[k+1])^2/3
  sigma[k+3,]=x
  sigma
}

#' Increasing concave function
#'
#' @param x value of the covariate
#' @param t value of t
#' @return value
#' @keywords CondCop
#' @export
incconcave=function(x,t){
  n=length(x)
  k=length(t)-2
  m=k+3
  sigma=matrix(1:(m*n)*0,nrow=m,ncol=n)
  for(j in 1:k){
    i1=x<=t[j]
    sigma[j,i1] = x[i1]-t[1]
    i2=x>t[j]&x<=t[j+1]
    sigma[j,i2] = t[j]-t[1]+((t[j+1]-t[j])^3-(t[j+1]-x[i2])^3)/3/(t[j+1]-t[j])/(t[j+2]-t[j]) +(x[i2]-t[j])*(t[j+2]-t[j+1])/(t[j+2]-t[j])
    i3=x>t[j+1]&x<=t[j+2]
    sigma[j,i3] = t[j]-t[1] + (t[j+1]-t[j])^2/3/(t[j+2]-t[j]) + (t[j+2]-t[j+1])*(t[j+1]-t[j])/(t[j+2]-t[j]) +((t[j+2]-t[j+1])^3-(t[j+2]-x[i3])^3)/3/(t[j+2]-t[j+1])/(t[j+2]-t[j])
    i4=x>=t[j+2]
    sigma[j,i4] = t[j]-t[1] + (t[j+1]-t[j])^2/3/(t[j+2]-t[j]) + (t[j+2]-t[j+1])*(t[j+1]-t[j])/(t[j+2]-t[j]) +(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
  }
  i1=x<=t[2]
  sigma[k+1,i1]=-(t[2]-x[i1])^3/3/(t[2]-t[1])^2
  i2=x>t[2]
  sigma[k+1,i2]=0
  i1=x<=t[k+1]
  sigma[k+2,i1]=x[i1]-t[1]
  i2=x>t[k+1]&x<=t[k+2]
  sigma[k+2,i2]=t[k+1]-t[1]+((t[k+2]-t[k+1])^2*(x[i2]-t[k+1])-(x[i2]-t[k+1])^3/3)/(t[k+2]-t[k+1])^2
  i3=x>t[k+2]
  sigma[k+2,i3]=t[k+1]-t[1]+((t[k+2]-t[k+1])^2*(t[k+2]-t[k+1])-(t[k+2]-t[k+1])^3/3)/(t[k+2]-t[k+1])^2
  sigma[k+3,]=x
  sigma
}

#' Decreasing convex function
#'
#' @param x value of the covariate
#' @param t value of t
#' @return value
#' @keywords CondCop
#' @export
decconvex=function(x,t){
  n=length(x)
  k=length(t)-2
  m=k+3
  sigma=matrix(1:(m*n)*0,nrow=m,ncol=n)
  for(j in 1:k){
    i1=x<=t[j]
    sigma[j,i1] = x[i1]-t[1]
    i2=x>t[j]&x<=t[j+1]
    sigma[j,i2] = t[j]-t[1]+((t[j+1]-t[j])^3-(t[j+1]-x[i2])^3)/3/(t[j+1]-t[j])/(t[j+2]-t[j]) +(x[i2]-t[j])*(t[j+2]-t[j+1])/(t[j+2]-t[j])
    i3=x>t[j+1]&x<=t[j+2]
    sigma[j,i3] = t[j]-t[1] + (t[j+1]-t[j])^2/3/(t[j+2]-t[j]) + (t[j+2]-t[j+1])*(t[j+1]-t[j])/(t[j+2]-t[j]) +((t[j+2]-t[j+1])^3-(t[j+2]-x[i3])^3)/3/(t[j+2]-t[j+1])/(t[j+2]-t[j])
    i4=x>=t[j+2]
    sigma[j,i4] = t[j]-t[1] + (t[j+1]-t[j])^2/3/(t[j+2]-t[j]) + (t[j+2]-t[j+1])*(t[j+1]-t[j])/(t[j+2]-t[j]) +(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
  }
  i1=x<=t[2]
  sigma[k+1,i1]=-(t[2]-x[i1])^3/3/(t[2]-t[1])^2
  i2=x>t[2]
  sigma[k+1,i2]=0
  i1=x<=t[k+1]
  sigma[k+2,i1]=x[i1]-t[1]
  i2=x>t[k+1]&x<=t[k+2]
  sigma[k+2,i2]=t[k+1]-t[1]+((t[k+2]-t[k+1])^2*(x[i2]-t[k+1])-(x[i2]-t[k+1])^3/3)/(t[k+2]-t[k+1])^2
  i3=x>t[k+2]
  sigma[k+2,i3]=t[k+1]-t[1]+((t[k+2]-t[k+1])^2*(t[k+2]-t[k+1])-(t[k+2]-t[k+1])^3/3)/(t[k+2]-t[k+1])^2
  sigma[k+3,]=x
  sigma=-sigma
  sigma
}

#' Decreasing concave function
#'
#' @param x value of the covariate
#' @param t value of t
#' @return value
#' @keywords CondCop
#' @export
decconcave=function(x,t){
  n=length(x)
  k=length(t)-2
  m=k+3
  sigma=matrix(1:m*n,nrow=m,ncol=n)
  for(j in 1:(k-1)){
    i1=x<=t[j]
    sigma[j,i1] = 0
    i2=x>t[j]&x<=t[j+1]
    sigma[j,i2] = (x[i2]-t[j])^3 / (t[j+2]-t[j]) / (t[j+1]-t[j])/3
    i3=x>t[j+1]&x<=t[j+2]
    sigma[j,i3] = x[i3]-t[j+1]-(x[i3]-t[j+2])^3/(t[j+2]-t[j])/(t[j+2]-t[j+1])/3+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
    i4=x>t[j+2]
    sigma[j,i4]=(x[i4]-t[j+1])+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
  }
  i1=x<=t[k]
  sigma[k,i1] = 0
  i2=x>t[k]&x<=t[k+1]
  sigma[k,i2] = (x[i2]-t[k])^3 / (t[k+2]-t[k]) / (t[k+1]-t[k])/3
  i3=x>t[k+1]
  sigma[k,i3] = x[i3]-t[k+1]-(x[i3]-t[k+2])^3/(t[k+2]-t[k])/(t[k+2]-t[k+1])/3+(t[k+1]-t[k])^2/3/(t[k+2]-t[k])-(t[k+2]-t[k+1])^2/3/(t[k+2]-t[k])
  i1=x<=t[2]
  sigma[k+1,i1]=x[i1]-t[1]+(t[2]-x[i1])^3/(t[2]-t[1])^2/3
  i2=x>t[2]
  sigma[k+1,i2]=x[i2]-t[1]
  i1=x<=t[k+1]
  sigma[k+2,i1]=0
  i2=x>t[k+1]
  sigma[k+2,i2]=(x[i2]-t[k+1])^3/(t[k+2]-t[k+1])^2/3
  sigma[k+3,]=x
  sigma=-sigma
  sigma
}


