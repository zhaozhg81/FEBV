source("../R/F_Functions.R")
library(REBayes)
library(Rmosek)


CVF <-function(cvid, X_train_AML, X_train_ALL, X_valid_AML, X_valid_ALL, bw="ucv")
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
  
  
  alpha <- 0.05

  bandwidth <- c(0.17,0.17,0.17)
  
  ljs   <- LJS(sSq,df)
  opt   <- OPT(sSq,df)
  smy   <- SMY(sSq,df)
  vsh   <- VSH(sSq,df)
  feb   <- fEB(sSq,df, bandwidth)
  fes   <- fES(sSq,df, bandwidth)
  reb   <- GVmix(sSq,array(df,G))$dy
  Feb   <- FEB(sSq,df)
##  feb <- Feb
##  fes <- Feb
  
  tauSq.sSq     <- tauSqfn(sSq,d,muhat,alpha)
  CI.sSq        <- CI.t(sSq,df,d,theta,sigmaSq,muhat,alpha)
  CI.ljs        <- CI.est(ljs,tauSq.sSq,d,theta,sigmaSq,muhat,alpha)
  CI.opt        <- CI.est(opt,tauSq.sSq,d,theta,sigmaSq,muhat,alpha)
  CI.smy        <- CI.est(smy,tauSq.sSq,d,theta,sigmaSq,muhat,alpha)
  CI.vsh        <- CI.est(vsh,tauSq.sSq,d,theta,sigmaSq,muhat,alpha)
  CI.feb        <- CI.est(feb,tauSq.sSq,d,theta,sigmaSq,muhat,alpha)
  CI.fes        <- CI.est(fes,tauSq.sSq,d,theta,sigmaSq,muhat,alpha)
  CI.reb        <- CI.est(reb,tauSq.sSq,d,theta,sigmaSq,muhat,alpha)
  CI.Feb        <- CI.est(Feb,tauSq.sSq,d,theta,sigmaSq,muhat,alpha)
  
  
  idvec=cbind(CI.sSq[,7],CI.ljs[,7],CI.opt[,7],CI.smy[,7],CI.vsh[,7],CI.feb[,7],CI.fes[,7],CI.reb[,7],CI.Feb[,7])
  idvec
}



concordance_leukemia <- function(numSim=500, PLOT=FALSE, bw="defaults")
{
  ## Reading Data
  ##---------------------------------------------------------------------------------------------------------------------------
  
  leukemiadata <- read.table("../realdata/Textdata/Leukemia_x.txt", quote="\"", comment.char="")
  class        <- read.delim("../realdata/Textdata/Leukemia_y.txt", header=FALSE) # 1: AML,   2: ALL
  leukemia     <- t(leukemiadata)
  
  n1  <- length(class[class==1]) # n1=47
  n2  <- length(class[class==2]) # n2=25
  n   <- n1+n2
  p   <- ncol(leukemia)  # number of Variables
  
  ## The function for calculating confidence intervals
  ##--------------------------------------------------------------------------------------------------------------------------

  if(PLOT==FALSE)
    {
      
      ## Simulation
      ##--------------------------------------------------------------------------------------------------------------------------
      
      Crate2 <- array(0, c(numSim,9))
      
      for(i in 1:numSim)
        {
          if( i%%1==0 ){
            print(date())
            print( paste( "number of simulation:", i, sep="") )
          }
          
          ## Divide data set into Training and Validation Set
          grp               <- class  # 1: AML, 2: ALL
          CV                <- c(rep(2,nrow(grp)))
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
          
          
          idvec1=CVF(1, X_train_AML, X_train_ALL, X_valid_AML, X_valid_ALL, bw)
          sSq1=idvec1[,1]
          ljs1=idvec1[,2]
          opt1=idvec1[,3]
          smy1=idvec1[,4]
          vsh1=idvec1[,5]
          feb1=idvec1[,6]
          fes1=idvec1[,7]
          reb1=idvec1[,8]
          Feb1=idvec1[,9]
          
          idvec2=CVF(2, X_train_AML, X_train_ALL, X_valid_AML, X_valid_ALL, bw)
          sSq2=idvec2[,1]
          ljs2=idvec2[,2]
          opt2=idvec2[,3]
          smy2=idvec2[,4]
          vsh2=idvec2[,5]
          feb2=idvec2[,6]
          fes2=idvec2[,7]
          reb2=idvec2[,8]
          Feb2=idvec2[,9]
          
          
          idsSq=ifelse((sSq1+sSq2)==1,1,0)
          idljs=ifelse((ljs1+ljs2)==1,1,0)
          idopt=ifelse((opt1+opt2)==1,1,0)
          idsmy=ifelse((smy1+smy2)==1,1,0)
          idvsh=ifelse((vsh1+vsh2)==1,1,0)
          idfeb=ifelse((feb1+feb2)==1,1,0)
          idfes=ifelse((fes1+fes2)==1,1,0)
          idreb=ifelse((reb1+reb2)==1,1,0)
          idFeb=ifelse((Feb1+Feb2)==1,1,0)
          
          Crate2[i,]=c(mean(idsSq),mean(idljs),mean(idopt),mean(idsmy),mean(idvsh),mean(idfeb),mean(idfes),mean(idreb),mean(idFeb))
        }
      
      ## Results
      ##------------------------------------------------------------------------------------------------------------------------------------
      
      mean   <-round(apply(Crate2,2,mean),3)
      std    <-round(apply(Crate2,2,sd  ),3)
      table2 <- cbind(mean,std)
      
      path        <- paste("./Data/")
      path1       <- paste("./Figure/")
      filenameR   <- paste(path,"Leukemia_Concordance_Rate.Rdata", sep="")
      filenameT   <- paste(path,"Leukemia_Concordance_Table.Rdata", sep="")
      save(Crate2, file = filenameR)
      save(table2, file = filenameT)
    }else{
      path        <- paste("./Data/")
      path1       <- paste("./Figure/")
      filenameR   <- paste(path,"Leukemia_Concordance_Rate.Rdata", sep="")
      filenameT   <- paste(path,"Leukemia_Concordance_Table.Rdata", sep="")
      load(filenameR)
      load(filenameT)
    }
  
  
  figurename1 <- paste(path1,"Leukemia_Concordance.pdf",sep="")
  pdf(figurename1,width=12,height=5)
  boxplot(Crate2[,1],Crate2[,2],Crate2[,3],Crate2[,4],Crate2[,5],Crate2[,6],Crate2[,7],Crate2[,8],Crate2[,9],main="Rate of discordant pairs : Leukemia data", ylim=c(0.05,0.4), names=c("S^2","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","REBayes","F-EBV"),col=c("cyan3","blue","purple","green4","gold4", "brown","orange","green","red"),cex.lab=1.2,cex.axis=1.2,cex.main=1.5)
  dev.off()
  
  path1     <-"/home/zhaozhg/Dropbox/Apps/Overleaf/On F-modelling based Empirical Bayes Estimation of Variances/figure/"
  figurename1 <- paste(path1,"Leukemia_Concordance.pdf",sep="")
  pdf(figurename1,width=12,height=5)
  boxplot(Crate2[,1],Crate2[,2],Crate2[,3],Crate2[,4],Crate2[,5],Crate2[,6],Crate2[,7],Crate2[,8],Crate2[,9],main="Rate of discordant pairs : Leukemia data", ylim=c(0.05,0.4), names=c("S^2","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","REBayes","F-EBV"),col=c("cyan3","blue","purple","green4","gold4", "brown","orange","green","red"),cex.lab=1.2,cex.axis=1.2,cex.main=1.5)
  dev.off()
  
  
  boxplot(Crate2[,1],Crate2[,2],Crate2[,3],Crate2[,4],Crate2[,5],Crate2[,6],Crate2[,7],Crate2[,8],Crate2[,9],main="Rate of discordant pairs : Leukemia data", ylim=c(0.05,0.4), names=c("S^2","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","REBayes","F-EBV"),col=c("cyan3","blue","purple","green4","gold4", "brown","orange","green","red"),cex.lab=1.2,cex.axis=1.2,cex.main=1.5)
  
  ##---------------------------------------------------------------------------------------------------------------------------------
}
