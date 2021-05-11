## Confidence Interval : EB Version
##
##
## This function calcualte the confidence interval based on F-EBV estimator of variances.
## 
## argument:
##      sSq: Sample variance
##        df: degrees of freedom
##        x: the observed sample mean
##      alpha: coverage probability
##
##
## output:
##    interval: 
##      $lb : lower bound
##      $ub : upper bound


CI_FEB <- function(sSq,df,x,alpha=0.05){
  
  p         <- length(sSq)
  estsigmaSq <- FEB(sSq, df)
  
  muhat = mean(x)
  
  zalpha  <- qnorm(1-alpha,0,1)
  
  tauSqhat<- max((zalpha^2+zalpha*sqrt(zalpha^2+2*sum(sSq)))/p, mean((x-muhat)^2)-mean(sSq))
  
  zalpha2   <- qnorm(1-alpha/2,0,1)
  esttheta  <- (esttauSq)/(estsigmaSq+esttauSq)*xi+(estsigmaSq)/(estsigmaSq+esttauSq)*muhat
  lb        <- esttheta - sqrt( zalpha2^2 - log( esttauSq/(esttauSq+estsigmaSq) ) )*sqrt((esttauSq*estsigmaSq)/(estsigmaSq+esttauSq))
  ub        <- esttheta + sqrt( zalpha2^2 - log( esttauSq/(esttauSq+estsigmaSq) ) )*sqrt((esttauSq*estsigmaSq)/(estsigmaSq+esttauSq))
  length    <- (ub-lb)/2                 # length
  interval  <- cbind(lb,ub)
  interval

  }