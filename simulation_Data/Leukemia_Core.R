library(latex2exp)
library(Matrix)
library(REBayes)
library(Rmosek)


source("../R/F_Functions.R")

set.seed(2)

if(PLOT==FALSE)
  {
    
    ##-- Setting Parameters ---------------------------------------------------------------------------------------------------------------------------
    
    
    ## leukemia Data
    leukemiadata <- read.table("./Textdata/Leukemia_x.txt", quote="\"", comment.char="")
    class        <- read.delim("./Textdata/Leukemia_y.txt", header=FALSE) # 1: AML,   2: ALL
    leukemia     <- leukemiadata
    
    ##-- Setting Parameters ----------------------------------------------------------------------------------------------------------------------------------------------
    
    number.gene       <- nrow(leukemia)   # Dimension
    nobs    <- ncol(leukemia)

    gene.ind <- sample( c(1:number.gene), 2000, replace=FALSE )

    G <- 2000
    
    ##-- Common Core ----------------------------------------------------------------------------------------------------------------------------------
    
    df      <- 2*(n-1)   # Degree of Freedom
    mu      <- 0     # Mean of prior distrubution
    alpha   <- 0.05
    
    Ncol    <- 13
    theta   <- array(0, G)
    sigmaSq <- array(0, G)
    bon     <- array(0, G)
    sSq     <- array(0, G)
    ljs     <- array(0, G)
    opt     <- array(0, G)
    smy     <- array(0, G)
    vsh     <- array(0, G)
    feb     <- array(0, G)
    fes     <- array(0, G)
    reb     <- array(0, G)
    Feb     <- array(0, G)
          
    CI.bon.sel.unknown       <- array(0, c(numSim, N,Ncol))
    CI.sSq.sel.unknown       <- array(0, c(numSim, N,Ncol))
    CI.ljs.sel.unknown       <- array(0, c(numSim, N,Ncol))
    CI.opt.sel.unknown       <- array(0, c(numSim, N,Ncol))
    CI.smy.sel.unknown       <- array(0, c(numSim, N,Ncol))
    CI.vsh.sel.unknown       <- array(0, c(numSim, N,Ncol))
    CI.feb.sel.unknown       <- array(0, c(numSim, N,Ncol))
    CI.fes.sel.unknown       <- array(0, c(numSim, N,Ncol))
    CI.reb.sel.unknown       <- array(0, c(numSim, N,Ncol))
    CI.Feb.sel.unknown       <- array(0, c(numSim, N,Ncol))
    
    
    tauSqhat                 <- array(0, numSim)
    tauSqhatFeb              <- array(0, numSim)
    
    Result                   <- array(0, c(N,41))
    Result.avg.unknown       <- array(0, c(N,41))
    
    tauSqs   <- seq(0,3,3/N)[2:(N+1)]
    
    
    ##-- Simulation ---------------------------------------------------------------------------------------------------------------------------
    for(i in  1:numSim){
      if( i%% 10==0 ){
        print(date())
        print( paste( "number of simulation:", i, sep="") )
      }
      rand0   <- sample(1:nobs, size=2*n, replace=FALSE)  
      rand1   <- rand0[1:n]
      rand2   <- rand0[(n+1):(2*n)]

      leuk1  <- leukemia[gene.ind,rand1]       # 2000 x n
      leuk2  <- leukemia[gene.ind,rand2]       # 2000 x n
      mean1  <- apply(leuk1,1,mean)
      mean2  <- apply(leuk2,1,mean)
      s1Sq   <- apply(leuk1,1,var)
      s2Sq   <- apply(leuk2,1,var)
      sSq    <- (s1Sq + s2Sq)/2

            
      bandwidth <- c(0.17,0.17,0.17)
    
      
      bon        <- sSq*2/n
      ljs        <- LJS(sSq,df)*2/n
      opt        <- OPT(sSq,df)*2/n
      smy        <- SMY(sSq,df)*2/n
      vsh        <- VSH(sSq,df)*2/n
      feb        <- fEB(sSq,df,bandwidth)*2/n
      fes        <- fES(sSq,df,bandwidth)*2/n
      reb        <- GVmix(sSq,array(df,G))$dy*2/n
      Feb        <- FEB(sSq,df)*2/n
      ## feb = Feb
      ## fes = Feb
      ## reb = Feb
      ## vsh = Feb

      for(t in 1:N)
        {
          theta   <- rnorm(G,0,sqrt(tauSqs[t])) 
          xi      <- mean1 - mean2 + theta
          muhat   <- mean(xi) 
          T.stat  <- xi/sqrt(sSq)
          max.ind <- order(abs(T.stat),decreasing=TRUE)[1:1]
          
          tauSqhat[i]<- tauSqfn(sSq*2/n,xi,muhat,alpha)
          
          ## # CI : tauSq unknown
          CI.bon.unknown   <- CI.bonferroni(sSq,df,xi,theta,sigmaSq,muhat,alpha)
          CI.sSq.unknown   <- CI.t(sSq,df,xi,theta,sigmaSq,muhat,alpha)
          CI.ljs.unknown   <- CI.est(ljs,tauSqhat[i],xi,theta,sigmaSq,muhat,alpha)
          CI.opt.unknown   <- CI.est(opt,tauSqhat[i],xi,theta,sigmaSq,muhat,alpha)
          CI.smy.unknown   <- CI.est(smy,tauSqhat[i],xi,theta,sigmaSq,muhat,alpha)
          CI.vsh.unknown   <- CI.est(vsh,tauSqhat[i],xi,theta,sigmaSq,muhat,alpha)
          CI.feb.unknown   <- CI.est(feb,tauSqhat[i],xi,theta,sigmaSq,muhat,alpha)
          CI.fes.unknown   <- CI.est(fes,tauSqhat[i],xi,theta,sigmaSq,muhat,alpha)
          CI.reb.unknown   <- CI.est(reb,tauSqhat[i],xi,theta,sigmaSq,muhat,alpha)
          CI.Feb.unknown   <- CI.est(Feb,tauSqhat[i],xi,theta,sigmaSq,muhat,alpha)
          
          ## CI : tauSq unknown - selected
          CI.bon.sel.unknown[i,t,]   <- CI.bon.unknown[max.ind,]
          CI.sSq.sel.unknown[i,t,]   <- CI.sSq.unknown[max.ind,]
          CI.ljs.sel.unknown[i,t,]   <- CI.ljs.unknown[max.ind,]
          CI.opt.sel.unknown[i,t,]   <- CI.opt.unknown[max.ind,]
          CI.smy.sel.unknown[i,t,]   <- CI.smy.unknown[max.ind,]
          CI.vsh.sel.unknown[i,t,]   <- CI.vsh.unknown[max.ind,]
          CI.feb.sel.unknown[i,t,]   <- CI.feb.unknown[max.ind,]
          CI.fes.sel.unknown[i,t,]   <- CI.fes.unknown[max.ind,]
          CI.reb.sel.unknown[i,t,]   <- CI.reb.unknown[max.ind,]
          CI.Feb.sel.unknown[i,t,]   <- CI.Feb.unknown[max.ind,]
                    
        }
    }
    for(t in 1:N)
      {
          bon.sel.unknown.Result=apply(CI.bon.sel.unknown[,t,7:10],2,mean);
          sSq.sel.unknown.Result=apply(CI.sSq.sel.unknown[,t,7:10],2,mean);
          ljs.sel.unknown.Result=apply(CI.ljs.sel.unknown[,t,7:10],2,mean);
          opt.sel.unknown.Result=apply(CI.opt.sel.unknown[,t,7:10],2,mean);
          smy.sel.unknown.Result=apply(CI.smy.sel.unknown[,t,7:10],2,mean);
          vsh.sel.unknown.Result=apply(CI.vsh.sel.unknown[,t,7:10],2,mean);
          feb.sel.unknown.Result=apply(CI.feb.sel.unknown[,t,7:10],2,mean);
          fes.sel.unknown.Result=apply(CI.fes.sel.unknown[,t,7:10],2,mean);
          reb.sel.unknown.Result=apply(CI.reb.sel.unknown[,t,7:10],2,mean);
          Feb.sel.unknown.Result=apply(CI.Feb.sel.unknown[,t,7:10],2,mean);
          
          
          Result[t,]=c(tauSqs[t],bon=bon.sel.unknown.Result,sSq=sSq.sel.unknown.Result,ljs=ljs.sel.unknown.Result,opt=opt.sel.unknown.Result,smy=smy.sel.unknown.Result,vsh=vsh.sel.unknown.Result,feb=feb.sel.unknown.Result,fes=fes.sel.unknown.Result,reb=reb.sel.unknown.Result, Feb=Feb.sel.unknown.Result)
        }

    ##----------------------------------------------------------------------------------------------------------------------------------------------
    ##-- Saving RData -------------------------------------------------------------------------------------------------------------------
    
    path      <- paste("./Data/")
    path1     <- paste("./Figure/")
    settname  <- paste("Leukemia_n_",n,sep="")
    filename  <- paste(path, settname, ".Rdata", sep="")
    save(Result, file = filename)    
}
if(PLOT==TRUE)
  {
    ##-- Saving RData -------------------------------------------------------------------------------------------------------------------
    
    path      <- paste("./Data/")
    path1     <- paste("./Figure/")
    settname  <- paste("Leukemia_n_",n,sep="")
    filename  <- paste(path, settname, ".Rdata", sep="")
    load(filename)
    source("SavingFigures_Luek.R")
  }
