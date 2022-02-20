library(usethis)
library(devtools)
library(ashr)
library(SQUAREM)
library(qvalue)
library(vashr)
library(limma)
library(matrixcalc)
library(kedd)
library(numDeriv)
library(Matrix)
library(REBayes)
library(statmod)
source("../R/F_EBV.R")

#-- begin : functions ----------------------------------------------------------------------------------------------------------------------------------------------

GenerateVar <- function(G, prior, a, b)
{
  temp       <- runif(G)
  if(prior==1)
    sigmaSq    <- 1/rgamma(G,abs(a),b)
  if(prior==2)
    sigmaSq <- rlnorm(G,a,b)
  if(prior==4)
    sigmaSq <- (temp<0.2)*rlnorm(G,a,b)+(temp>0.2)*(temp<0.6)*rlnorm(G,0,sqrt(2*log(2))) + (temp>0.6)*rlnorm(G, 0, sqrt(4*log(2)))
  if(prior==3)
    sigmaSq <- (temp<0.2)*1/rgamma(G,a,b)+(temp>0.2)*(temp<0.6)*1/rgamma(G,8,6) + (temp>0.6)*1/rgamma(G,9,19)
  if(prior==5)
    sigmaSq <- rinvgauss(G, 1/a, 1)
  if(prior==6)
    sigmaSq <- (temp<0.4) * rinvgauss(G, 1/a, 1) + (temp>0.4) * rinvgauss(G, a, a^4)
  if(prior==7)
    sigmaSq <- (temp<0.4) * rep(1/a, G) + (temp>0.4) * rep(a, G)

    
  sigmaSq
}

LJST   <- function(sSq,df){
  p    <- length(sSq)
  sSqp <- prod(sSq^(1/p))
  m    <- digamma(df/2)-log(df/2) # mean(log(chi2/df))
  V    <- trigamma(df/2)          # var(log(chi2/df))
  B    <- exp(-m)
  a0   <- max(0,1-(p-3)*V/sum((log(sSq)-mean(log(sSq)))^2))
  ljs  <- B*(sSqp^(1-a0))*(sSq^a0)
  c(B,a0)
}

OPTT       <- function(sSq,df){
  y        <- c(seq(0,1,0.01))
  G        <- length(sSq)
  sSqpool  <- exp(sum(log(sSq))/G)
  
  sortsSq  <- sort(sSq)
  sSq1     <- sortsSq[((G*0.1)+1):G]
  G1       <- length(sSq1)
  sSqpool1 <- exp(sum(log(sSq1))/G)
  q        <- sSqpool1/sSq1
  r1       <- matrix(rep(q,each=length(y))^(y)  ,length(y),length(q))
  r2       <- matrix(rep(q,each=length(y))^(2*y),length(y),length(q))
  rm1      <- apply(r1,1,mean)
  rm2      <- apply(r2,1,mean)
  risk     <- rm2*(hf(df,1,G1)^(2*y))*(hf(df,1,1)^(2*(1-y)))/(hf(df,2*y/G1,1)^(G1-1))/(hf(df,2*(1-y+y/G1),1))-(2*rm1)*(hf(df,1,G1)^(y))*(hf(df,1,1)^(1-y))/(hf(df,y/G1,1)^(G1-1))/(hf(df,(1-y+y/G1),1))+1
  ahat     <- y[which.min(risk)]
  opt_chqbc<- ((hf(df,1,G)*sSqpool)^ahat)*(hf(df,1,1)*sSq^(1-ahat))
  
  nsSqpool  <- exp(sum(log(opt_chqbc))/G)
  sortsSq0  <- sort(opt_chqbc)
  sSq0      <- sortsSq0[((G*0.1)+1):G]
  G0        <- length(sSq0)
  nsSqpool0 <- exp(sum(log(sSq0))/G)
  nq        <- sSq0/nsSqpool0
  nr1       <- matrix(rep(nq,each=length(y))^(y)  ,length(y),length(nq))
  nr2       <- matrix(rep(nq,each=length(y))^(2*y),length(y),length(nq))
  nrm1      <- apply(nr1,1,mean)
  nrm2      <- apply(nr2,1,mean)
  nrisk     <- nrm2*hf(df,1,G0)^(2*y)*hf(df,1,1)^(2*(1-y))/hf(df,2*y/G0,1)^(G0-1)/hf(df,2*(1-y+y/G0),1)-2*nrm1*hf(df,1,G0)^(y)*hf(df,1,1)^(1-y)/hf(df,y/G0,1)^(G0-1)/hf(df,(1-y+y/G0),1)+1
  nahat     <- y[which.min(nrisk)]
  opt.chqbc <- ((hf(df,1,G)*nsSqpool)^nahat)*(hf(df,1,1)*opt_chqbc^(1-nahat))
  c(sSqpool,ahat,nsSqpool,nahat)
}




## Koencker method

## Exponential Lindley-James-Stein Estimator
LJS    <- function(sSq,df){
  p    <- length(sSq)
  sSqp <- prod(sSq^(1/p))
  m    <- digamma(df/2)-log(df/2) # mean(log(chi2/df))
  V    <- trigamma(df/2)          # var(log(chi2/df))
  B    <- exp(-m)
  a0   <- max(0,1-(p-3)*V/sum((log(sSq)-mean(log(sSq)))^2))
  ljs  <- B*(sSqp^(1-a0))*(sSq^a0)
  ljs
}

## f-EBV
fEB     <- function(sSq,df, bandwidth){
  f0    <-dkde(sSq,deriv.order=0, h=bandwidth[1])
  f1    <-dkde(sSq,deriv.order=1, h=bandwidth[2])
  fval0 <-approx(f0$eval.points, f0$est.fx, sSq)
  fval1 <-approx(f1$eval.points, f1$est.fx, sSq)
  feb   <-((df-2)/(df*sSq)-2/df*(fval1$y/max(fval0$y,0.00001)))^{-1}
  feb   <-ifelse(feb<0,sSq,feb)
  feb
}

fES     <- function(sSq,df, bandwidth){
  f0    <-dkde(sSq,deriv.order=0, h=bandwidth[1])
  f1    <-dkde(sSq,deriv.order=1, h=bandwidth[2])
  f2    <-dkde(sSq,deriv.order=2, h=bandwidth[3])
  fval0 <-approx(f0$eval.points, f0$est.fx, sSq)
  fval1 <-approx(f1$eval.points, f1$est.fx, sSq)
  fval2 <-approx(f2$eval.points, f2$est.fx, sSq)
  febvs <-(df*(df-2)*sSq*fval0$y-2*df*sSq^2*fval1$y)/(4*sSq^2*fval2$y-4*(df-2)*sSq*fval1$y+df*(df-2)*fval0$y)
  febvs <-ifelse(febvs<0,sSq,febvs)
  febvs
}


## h function
hf <- function(df,t,p){
  h <-(df/2)^t*(gamma(df/2)/gamma(df/2+t/p))^p
  h
}

## Optimal Estimator (Tong)
OPT        <- function(sSq,df){
  y        <- c(seq(0,1,0.01))
  G        <- length(sSq)
  sSqpool  <- exp(sum(log(sSq))/G)
  
  sortsSq  <- sort(sSq)
  sSq1     <- sortsSq[((G*0.1)+1):G]
  G1       <- length(sSq1)
  sSqpool1 <- exp(sum(log(sSq1))/G)
  q        <- sSqpool1/sSq1
  r1       <- matrix(rep(q,each=length(y))^(y)  ,length(y),length(q))
  r2       <- matrix(rep(q,each=length(y))^(2*y),length(y),length(q))
  rm1      <- apply(r1,1,mean)
  rm2      <- apply(r2,1,mean)
  risk     <- rm2*(hf(df,1,G1)^(2*y))*(hf(df,1,1)^(2*(1-y)))/(hf(df,2*y/G1,1)^(G1-1))/(hf(df,2*(1-y+y/G1),1))-(2*rm1)*(hf(df,1,G1)^(y))*(hf(df,1,1)^(1-y))/(hf(df,y/G1,1)^(G1-1))/(hf(df,(1-y+y/G1),1))+1
  ahat     <- y[which.min(risk)]
  opt_chqbc<- ((hf(df,1,G)*sSqpool)^ahat)*(hf(df,1,1)*sSq^(1-ahat))
  
  nsSqpool  <- exp(sum(log(opt_chqbc))/G)
  sortsSq0  <- sort(opt_chqbc)
  sSq0      <- sortsSq0[((G*0.1)+1):G]
  G0        <- length(sSq0)
  nsSqpool0 <- exp(sum(log(sSq0))/G)
  nq        <- sSq0/nsSqpool0
  nr1       <- matrix(rep(nq,each=length(y))^(y)  ,length(y),length(nq))
  nr2       <- matrix(rep(nq,each=length(y))^(2*y),length(y),length(nq))
  nrm1      <- apply(nr1,1,mean)
  nrm2      <- apply(nr2,1,mean)
  nrisk     <- nrm2*hf(df,1,G0)^(2*y)*hf(df,1,1)^(2*(1-y))/hf(df,2*y/G0,1)^(G0-1)/hf(df,2*(1-y+y/G0),1)-2*nrm1*hf(df,1,G0)^(y)*hf(df,1,1)^(1-y)/hf(df,y/G0,1)^(G0-1)/hf(df,(1-y+y/G0),1)+1
  nahat     <- y[which.min(nrisk)]
  opt.chqbc <- ((hf(df,1,G)*nsSqpool)^nahat)*(hf(df,1,1)*opt_chqbc^(1-nahat))
  opt.chqbc
}

## Smyth Estimator 
SMY    <- function(sSq,df){ 
  smy    <-squeezeVar(sSq, df)$var.post
  smy
}

modified.SMY    <- function(sSq,df){
  #smy    <-squeezeVar(sSq, df)$var.post
  #smy
  
  res    <-squeezeVar(sSq, df)
  a0 = res$df.prior/2 # prior parameters
  b0 = (res$df.prior/2)*res$var.prior
  a1 = a0 + df/2
  b1 = b0 + df*sSq/2
  smy <- b1/(a1-2) # instead of usual estimate, b1/a1, which is inverse of posterior mean of precision
  smy
}

VSH  <- function(sSq, df){      
  vsh <- (vash(sqrt(sSq),df)$sd.post)^2
  vsh
}  

modified.VSH <- function(sSq, df){
  G <- length(sSq)
  vsh <- vash( sqrt(sSq),df )
  a0 = vsh$fitted.g$alpha
  b0 = vsh$fitted.g$beta
  pi0 = matrix( rep( vsh$fitted.g$pi, G), c(length(a0), G) )
  
  a1 = matrix( rep( a0 + df/2, G), c(length(a0), G) )
  b1 = matrix( rep(b0,G), c(length(b0), G ) ) + df* t(matrix( rep(sSq,length(b0)), c(G,length(b0) ) ))/2

  vsh <- apply( pi0 * b1/(a1-2), 2, sum )
  vsh
}


## Estimator for tauSq
tauSqfn   <- function(estsigmaSq,xi,muhat,alpha=0.05,trunc=TRUE){
  G       <- length(estsigmaSq)
  zalpha  <- qnorm(1-alpha,0,1)
  if( trunc==TRUE){
    tauSqhat<- max((zalpha^2+zalpha*sqrt(zalpha^2+2*sum(estsigmaSq)))/G, mean((xi-muhat)^2)-mean(estsigmaSq))
  }else{
    tauSqhat<-  mean((xi-muhat)^2)-mean(estsigmaSq)
  }
  tauSqhat
}

## Confidence Interval : EB Version
CI.est <- function(estsigmaSq,esttauSq,xi,theta,sigmaSq,muhat,alpha){
  p         <- length(estsigmaSq)
  zalpha2   <- qnorm(1-alpha/2,0,1)
  esttheta  <- (esttauSq)/(estsigmaSq+esttauSq)*xi+(estsigmaSq)/(estsigmaSq+esttauSq)*muhat
  lb        <- esttheta - sqrt( zalpha2^2 - log( esttauSq/(esttauSq+estsigmaSq) ) )*sqrt((esttauSq*estsigmaSq)/(estsigmaSq+esttauSq))
  ub        <- esttheta + sqrt( zalpha2^2 - log( esttauSq/(esttauSq+estsigmaSq) ) )*sqrt((esttauSq*estsigmaSq)/(estsigmaSq+esttauSq))
  idx       <- (lb<=theta)*(theta<=ub)   # coverage
  length    <- (ub-lb)/2                 # length
  riskmean  <- (esttheta-theta)^2        # risk of mean
  risksigma <- (sigmaSq/estsigmaSq-1)^2  # risk of sigma
  interval  <- cbind(theta,esttheta,sigmaSq,estsigmaSq,lb,ub,idx,length,risksigma,riskmean,esttauSq,xi,muhat)
  interval
}

## Confidence Interval : for bonferroni
CI.bonferroni <- function(estsigmaSq,df,xi,theta,sigmaSq,muhat,alpha)
  {
    p        <- length(estsigmaSq)
    esttheta <- xi
    lb       <- esttheta - qt(1-alpha/(2*p),df)*sqrt(estsigmaSq)
    ub       <- esttheta + qt(1-alpha/(2*p),df)*sqrt(estsigmaSq)
    idx      <- (lb<=theta)*(theta<=ub)   # coverage
    length   <- (ub-lb)/2                 #length
    riskmean <- (esttheta-theta)^2        # risk of mean
    risksigma<- (sigmaSq/estsigmaSq-1)^2  # risk of sigma
    interval <- cbind(theta,esttheta,sigmaSq,estsigmaSq,lb,ub,idx,length,risksigma,riskmean,df,xi,muhat)
    interval
  }

## Confidence Interval: t
CI.t <- function(estsigmaSq,df,xi,theta,sigmaSq,muhat,alpha)
  {
    esttheta  <- xi
    lb        <- esttheta - qt(1-alpha/2,df)*sqrt(estsigmaSq)
    ub        <- esttheta + qt(1-alpha/2,df)*sqrt(estsigmaSq)
    idx       <- (lb<=theta)*(theta<=ub)   # coverage
    length    <- (ub-lb)/2                 # length
    riskmean  <- (esttheta-theta)^2        # risk of mean
    risksigma <- (sigmaSq/estsigmaSq-1)^2  # risk of sigma
    interval  <- cbind(theta,esttheta,sigmaSq,estsigmaSq,lb,ub,idx,length,risksigma,riskmean,df,xi,muhat)
    interval
  }

changenames <- function(mtx){
  colnames(mtx,do.NULL=FALSE)
  colnames(mtx)<-c("theta","esttheta","sigmaSq","estsigmaSq","lb","ub","idx","length","risksigma","riskmean","esttauSq","xi","muhat")
  mtx
}

newton.raphson <- function(f, a, b, tol = 1e-5, n = 1000) {
  require(numDeriv) # Package for computing f'(x)
  
  x0 <- a # Set start value to supplied lower bound
  k  <- n # Initialize for iteration results
  
                                        # Check the upper and lower bounds to see if approximations result in 0
  fa <- f(a)
  if (fa == 0.0) {
    return(a)
  }
  
  fb <- f(b)
  if (fb == 0.0) {
    return(b)
  }
  
  for (i in 1:n) {
    dx   <- genD(func = f, x = x0)$D[1] # First-order derivative f'(x0)
    x1   <- x0 - (f(x0) / dx) # Calculate next value x1
    k[i] <- x1 # Store x1
    ## Once the difference between x0 and x1 becomes sufficiently small, output the results.
    if (abs(x1 - x0) < tol) {
      root.approx <- tail(k, n=1)
      res <- list('root approximation' = root.approx, 'iterations' = k)
      return(res)
    }
    ## If Newton-Raphson has not yet reached convergence set x1 as x0 and continue
    x0 <- x1
  }
  print('Too many iterations in method')
}

##-- end : functions ----------------------------------------------------------------------------------------------------------------------------------------------
