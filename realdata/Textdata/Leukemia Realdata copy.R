
source("../F_Functions.R")

# Reading Data
#---------------------------------------------------------------------------------------------------------------------------

leukemiadata <- read.table("./Leukemia_x.txt", quote="\"", comment.char="")
class        <- read.delim("./Leukemia_y.txt", header=FALSE) # 1: AML,   2: ALL
leukemia     <- t(leukemiadata)

n1  <- length(class[class==1]) # n1=47
n2  <- length(class[class==2]) # n2=25
n   <- n1+n2
p   <- ncol(leukemia)  # number of Variables

# The function for calculating confidence intervals
#--------------------------------------------------------------------------------------------------------------------------
CVF <-function(cvid)
{

  if (cvid==1)
  {
    con <-X_train_AML
    exp <-X_train_ALL
  }
  else if (cvid==2)
  {
    con <-X_valid_AML
    exp <-X_valid_ALL
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

  sigmaSq.ljs   <- array(0, G)
  sigmaSq.opt   <- array(0, G)
  sigmaSq.feb   <- array(0, G)
  sigmaSq.Feb   <- array(0, G)

  Ncol          <- 13
  CI.sSq        <- array(0, c(G,Ncol))
  CI.ljs        <- array(0, c(G,Ncol))
  CI.opt        <- array(0, c(G,Ncol))
  CI.feb        <- array(0, c(G,Ncol))
  CI.Feb        <- array(0, c(G,Ncol))

  alpha         <- 0.05
  sigmaSq.ljs   <- LJS(sSq,df)
  sigmaSq.opt   <- OPT(sSq,df)
  sigmaSq.feb   <- fEB(sSq,df)
  sigmaSq.Feb   <- FEB(sSq,df)

  tauSq.sSq     <- tauSqfn(sSq,d,muhat,alpha)
  CI.sSq        <- CI.t(sSq,df,d,theta,sigmaSq,muhat,alpha)
  CI.ljs        <- CI.est(sigmaSq.ljs,tauSq.sSq,d,theta,sigmaSq,muhat,1,alpha)
  CI.opt        <- CI.est(sigmaSq.opt,tauSq.sSq,d,theta,sigmaSq,muhat,1,alpha)
  CI.feb        <- CI.est(sigmaSq.feb,tauSq.sSq,d,theta,sigmaSq,muhat,1,alpha)
  CI.Feb        <- CI.est(sigmaSq.Feb,tauSq.sSq,d,theta,sigmaSq,muhat,1,alpha)

  idvec=cbind(CI.sSq[,7],CI.ljs[,7],CI.opt[,7],CI.feb[,7],CI.Feb[,7])
  idvec
}


# Simulation
#--------------------------------------------------------------------------------------------------------------------------
numSim=1000
Crate1 <- array(0, c(numSim,5))

for(i in 1:numSim)
{
  if( i%%10==0 ){
    print(date())
    print( paste( "number of simulation:", i, sep="") )
  }

  ## Divide data set into Training and Validation Set
  grp                <- class  # 1: AML, 2: ALL
  CV                 <- c(rep(2,nrow(grp)))
  train_AML_ind     <- sample(seq(length(grp[grp==1])+1,nrow(grp),1), size=length(grp[grp==1])/2) # Sampling Traing set
  train_ALL_ind     <- sample(seq(1,length(grp[grp==2]),1), size=length(grp[grp==2])/2)           # Sampling Traing set
  CV[train_AML_ind] <- 1
  CV[train_ALL_ind] <- 1
  Xdata=cbind(CV, grp, leukemia)

  ## CV - 1: Traing, 2: Validation;  grp - 1: AML, 2: ALL
  X_train_AML  <- Xdata[CV==1 & grp==1,3:(p+2)]
  X_train_ALL  <- Xdata[CV==1 & grp==2,3:(p+2)]
  X_valid_AML  <- Xdata[CV==2 & grp==1,3:(p+2)]
  X_valid_ALL  <- Xdata[CV==2 & grp==2,3:(p+2)]

  idvec1=CVF(1)
  sSq1=idvec1[,1]
  ljs1=idvec1[,2]
  opt1=idvec1[,3]
  feb1=idvec1[,4]
  Feb1=idvec1[,5]

  idvec2=CVF(2)
  sSq2=idvec2[,1]
  ljs2=idvec2[,2]
  opt2=idvec2[,3]
  feb2=idvec2[,4]
  Feb2=idvec2[,5]

  idsSq=ifelse((sSq1+sSq2)==1,1,0)
  idljs=ifelse((ljs1+ljs2)==1,1,0)
  idopt=ifelse((opt1+opt2)==1,1,0)
  idfeb=ifelse((feb1+feb2)==1,1,0)
  idFeb=ifelse((Feb1+Feb2)==1,1,0)

  Crate1[i,]=c(mean(idsSq),mean(idljs),mean(idopt),mean(idfeb),mean(idFeb))
}

# Results
#------------------------------------------------------------------------------------------------------------------------------------

mean<-round(apply(Crate1,2,mean),3)
std <-round(apply(Crate1,2,sd  ),3)
cbind(mean,std)
boxplot(Crate1[,1],Crate1[,2],Crate1[,3],Crate1[,4],Crate1[,5],main="Rate of discordant pairs : Leukemia data", ylim=c(0.05,0.4), names=c("S^2","ELJS","TW","fEBV","FEBV"),col=c("Green","Blue","Purple","Brown","Red"),cex.lab=1.2,cex.axis=1.2,cex.main=1.5)
#---------------------------------------------------------------------------------------------------------------------------------
