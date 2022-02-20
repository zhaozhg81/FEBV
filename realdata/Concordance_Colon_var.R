  source("../R/F_Functions.R")
library(REBayes)
library(Rmosek)

CVF <-function(cvid, X_train_tumor, X_train_normal, X_valid_tumor, X_valid_normal, bw="ucv")
{
  
  if (cvid==1)
  {
    con <-X_train_tumor
    exp <-X_train_normal
  }
  else if (cvid==2)
  {
    con <-X_valid_tumor
    exp <-X_valid_normal
  }
  
  m1    <- nrow(con)
  m2    <- nrow(exp)
  m     <- m1+m2
  c0    <- sqrt(1/m1+1/m2)
  G     <- ncol(con)  # Dimension
  df    <- m-2        # Digree of freedom
  
  x1    <- apply(con, 2, mean)
  x2    <- apply(exp, 2, mean)
  d     <- x2-x1  # mean difference #
  
  s1    <- apply(con, 2, var)
  s2    <- apply(exp, 2, var)
  sp    <- ((m1-1)*s1+(m2-1)*s2)/(m1+m2-2)
  sSq   <- (c0^2)*sp
  muhat <- mean(d)
  
  sigmaSq       <- array(0, G)
  theta         <- array(0, G)
  
  
  alpha <- 0.05
  bandwidth <- array(0.6, 3)
  if( bw == 'ucv' )
    {
      bandwidth[1] <- h.ucv(sSq, deriv.order=0)$h
      bandwidth[2] <- h.ucv(sSq, deriv.order=1)$h
      bandwidth[3] <- h.ucv(sSq, deriv.order=2)$h                
    }
  if( bw=='defaults')
    {
      bandwidth=c(0.6,0.6,0.6)
    }
  
  ljs   <- LJS(sSq,df)
  opt   <- OPT(sSq,df)
  
  modified.smy   <- modified.SMY(sSq,df)
  smy   <- SMY(sSq,df)
  vsh   <- modified.VSH(sSq,df)
  #feb   <- fEB(sSq,df, bandwidth)
  #fes   <- fES(sSq,df, bandwidth)
  #reb   <- GVmix(sSq,array(df,G))$dy
  Feb   <- FEB(sSq,df)

  feb <- Feb
  fes <- Feb
  reb <- Feb

  tauSq.sSq     <- tauSqfn(sSq,d,muhat,alpha)
  
  idvec=list(sSq=sSq, ljs=ljs, opt=opt, smy=smy, modified.smy=modified.smy, vsh=vsh, feb=feb, fes=fes, reb=reb, Feb=Feb)
  idvec
}



##concordance_colon <- function(numSim=500, PLOT=FALSE, bw="ucv")
##  {

numSim <- 500
bw <- "defaults"

    ## Reading Data
    ##---------------------------------------------------------------------------------------------------------------------------
    coloncancer <- read.table("../realdata/Textdata/colon_x.txt", quote="\"", comment.char="")
    tissue      <- read.delim("../realdata/Textdata/colon_y.txt", header=FALSE) # 1: tumor,   2: normal
    colon       <- t(coloncancer)
    
    n1  <- length(tissue[tissue==1]) # n1=40
    n2  <- length(tissue[tissue==2]) # n2=22
    n   <- n1+n2
    p   <- ncol(colon)  # number of Variables
    

    
    ## Simulation
    ##--------------------------------------------------------------------------------------------------------------------------
    
    est.error <- array(0, c(numSim, 10) )
    
    for(numsim in 1:numSim)
      {
        if( numsim%%1==0 ){
          print(date())
          print( paste( "number of simulation:", numsim, sep="") )
        }
        
        ## Divide data set into Training and Validation Set
        grp                <- tissue  # 1: tumor, 2: normal
        CV                 <- c(rep(2,nrow(grp)))
        train_tumr_ind     <- sample(seq(length(grp[grp==1])+1,nrow(grp),1), size=length(grp[grp==1])/2) # Sampling Traing set
        train_norm_ind     <- sample(seq(1,length(grp[grp==2]),1), size=length(grp[grp==2])/2)           # Sampling Traing set
        CV[train_norm_ind] <- 1
        CV[train_tumr_ind] <- 1
        Xdata=cbind(CV, grp, colon)
            ## CV - 1: Traing, 2: Validation;  grp - 1: tumor, 2: normal
        X_train_tumor  <- Xdata[CV==1 & grp==1,3:(p+2)]   
        X_train_normal <- Xdata[CV==1 & grp==2,3:(p+2)]
        X_valid_tumor  <- Xdata[CV==2 & grp==1,3:(p+2)]
        X_valid_normal <- Xdata[CV==2 & grp==2,3:(p+2)]
        
        
        var1.est = CVF(1, X_train_tumor, X_train_normal, X_valid_tumor, X_valid_normal, bw)                 
        var2.est = CVF(2, X_train_tumor, X_train_normal, X_valid_tumor, X_valid_normal, bw)        

        ORD <- order(var1.est$sSq)[1:50]

        est.error[numsim,] <- c( mean( (var2.est$sSq[ORD]/var1.est$sSq[ORD]-1)^2 ), mean( (var2.est$sSq[ORD]/var1.est$ljs[ORD]-1)^2 ), mean( (var2.est$sSq[ORD]/var1.est$opt[ORD]-1)^2 ),
                                mean( (var2.est$sSq[ORD]/var1.est$smy[ORD]-1)^2 ), mean( (var2.est$sSq[ORD]/var1.est$modified.smy[ORD]-1)^2 ),
                                mean( (var2.est$sSq[ORD]/var1.est$vsh[ORD]-1)^2 ), mean( (var2.est$sSq[ORD]/var1.est$feb[ORD]-1)^2 ),
                              mean( (var2.est$sSq[ORD]/var1.est$fes[ORD]-1)^2 ), mean( (var2.est$sSq[ORD]/var1.est$reb[ORD]-1)^2 ), mean( (var2.est$sSq[ORD]/var1.est$Feb[ORD]-1)^2 )
                              )
      }
##  }
