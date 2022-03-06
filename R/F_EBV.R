## This is the main routine to compute the F-modeling based empirical Bayes estimator of the variances (F-EBV)

FEB       <- function(sSq,df, K=10){
  ## The function take three arguments
  ## 1. sSq, the sample variances;
  ## 2. df, the degree of freedom. In this function, the degree of freedom is a constant across all the sample variances;
  ## 3. K, the parameters corresponding to the top K largest sample variances are estimator by the corresponding sample variance. By default, this value is taken as 10.

  ## Output:
  ## sigmaSq.hat, the FEBV estimator of the variances.
  
  p        <- length(sSq)
  ORD      <- order(sSq)
  sSq.sort <- sSq[ORD]
  sSqlag1  <- sSq.sort[2:p]
  sSqlag   <- c(sSqlag1,sSq.sort[p])
  hvec     <- sSq.sort/sSqlag
  hvec_num <- hvec^(df/2-2)
  hvec_den <- hvec^(df/2-1)
  rvec_num <- array(1,p)
  rvec_den <- array(1,p)
  
  for(i in 2:(p-1)){
    rvec_num[(p-i)]<- 1+hvec_num[(p-i+1)]*rvec_num[(p-i+1)]
    rvec_den[(p-i)]<- 1+hvec_den[(p-i+1)]*rvec_den[(p-i+1)]
  }
  rvec             <- rvec_num/rvec_den
  temp             <- df/2*(sSqlag*rvec-sSq.sort)
  temp[p]          <- NaN
  temp[(p-K+1):p]      <- sSq.sort[ (p-K+1):p ]
  sigmaSq.hat      <- temp
  sigmaSq.hat[ORD] <- temp
  sigmaSq.hat
}
