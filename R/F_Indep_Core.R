library(latex2exp)
library(REBayes)
library(Rmosek)

CI_Ind <- function(G, numSim, N, n, prior, a, b, PLOT=FALSE, bw= 'ucv' )
{

  
  if(PLOT==FALSE)
    {
      
      ## ##-- Setting Parameters ----------------------------------------------------------------------------------------------------------------------------------------------
      
      ## G       <- 1000   # Dimension
      ## numSim  <- 500   # The number of simulation
      ## N       <- 20     # The Number of ticks for X-axis
      ## n       <- 6      # Sample size  (df=n-1)

      ##-- Indep Core ------------------------------------------------------------------------------------------------------------------------------------------------------
      
      df      <- n-1   # Degree of Freedom
      mu      <- 0     # Mean of prior distrubution
      alpha   <- 0.05
      
      GAM     <-0
      INVGAM  <-0
      LNORM   <-0
      MIXIG   <-0
      MIXGAM  <-0

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

 
      bon.sel.unknown.Result   <- array(0, c(N, 4))
      sSq.sel.unknown.Result   <- array(0, c(N, 4))
      ljs.sel.unknown.Result   <- array(0, c(N, 4))
      opt.sel.unknown.Result   <- array(0, c(N, 4))
      smy.sel.unknown.Result   <- array(0, c(N, 4))
      vsh.sel.unknown.Result   <- array(0, c(N, 4))
      feb.sel.unknown.Result   <- array(0, c(N, 4))
      fes.sel.unknown.Result   <- array(0, c(N, 4))
      reb.sel.unknown.Result   <- array(0, c(N, 4))
      Feb.sel.unknown.Result   <- array(0, c(N, 4))
      
      bonR                     <- array(0, c(N, 4))
      sSqR                     <- array(0, c(N, 4))
      ljsR                     <- array(0, c(N, 4))
      optR                     <- array(0, c(N, 4))
      smyR                     <- array(0, c(N, 4))
      vshR                     <- array(0, c(N, 4))
      febR                     <- array(0, c(N, 4))
      fesR                     <- array(0, c(N, 4))
      rebR                     <- array(0, c(N, 4))
      FebR                     <- array(0, c(N, 4))
      
      Result           <- array(0, c(N,41))
      
      
      M=seq(0,1,0.999/N)
      M=M[2:(N+1)]
      
      
      ##-- Simulation ---------------------------------------------------------------------------------------------------------------------------------------------------------
      
      for(i in  1:numSim){
        if( i%% 10==0 )
          {
            print(date())
            print( paste( "number of simulation:", i, sep="") )
        }
        
        sigmaSq <- GenerateVar(G, prior, a, b)
        
        
        sSq        <- sigmaSq*rchisq(G,df)/df              
        bon        <- sSq
        sSq        <- sSq
        
        if( i==1 )
          {
            if( bw=="ucv")
              {
                bandwidth <- array(0, 3)
                bandwidth[1] <- h.ucv(sSq, deriv.order=0)$h
                bandwidth[2] <- h.ucv(sSq, deriv.order=1)$h
                bandwidth[3] <- h.ucv(sSq, deriv.order=2)$h
              }else{
                bandwidth <- array(0.2, 3)
              }
          }
        
        ljs        <- LJS(sSq,df)
        opt        <- OPT(sSq,df)
        smy        <- SMY(sSq,df)
        vsh        <- VSH(sSq,df)
        feb        <- fEB(sSq,df,bandwidth)
        fes        <- fES(sSq,df,bandwidth)
        reb        <- GVmix(sSq, array(df,G) )$dy
        Feb        <- FEB(sSq,df)
        
        ## tauSqhat.sSq[i]<- tauSqfn(sSq,xi,muhat,alpha)
        ## tauSqhat[i] <- tauSqfn(sSq,xi,muhat,alpha)
        ## tauSqhat.ljs[i]<- tauSqfn(ljs,xi,muhat,alpha)
        ## tauSqhat.opt[i]<- tauSqfn(opt,xi,muhat,alpha)
        ## tauSqhat.smy[i]<- tauSqfn(smy,xi,muhat,alpha)
        ## tauSqhat.vsh[i]<- tauSqfn(vsh,xi,muhat,alpha)
        ## tauSqhat.feb[i]<- tauSqfn(feb,xi,muhat,alpha)
        ## tauSqhat.fes[i]<- tauSqfn(fes,xi,muhat,alpha)
          ## tauSqhat.reb[i]<- tauSqfn(reb,xi,muhat,alpha)
        ## tauSqhat.Feb[i]<- tauSqfn(Feb,xi,muhat,alpha)
        
        
        for(t in 1:N)
          {
            
            tauSq      <-  M[t]/(1-M[t])
            theta      <- rnorm(G,mu,sqrt(tauSq))
            xi         <- rnorm(G,theta,sqrt(sigmaSq))
            muhat      <- mean(xi) #sum((xi/sSq))/sum(1/sSq)
            T.stat     <- xi/sqrt(sSq)
            max.ind    <- order(abs(T.stat),decreasing=TRUE)[(G/G):(G/G)]
            
            tauSqhat <- tauSqfn(sSq,xi,muhat,alpha)
            tauSqhat = tauSq 
            tauSqhat.sSq <- tauSqhat
            tauSqhat.ljs <- tauSqhat
            tauSqhat.opt <- tauSqhat
            tauSqhat.smy <- tauSqhat
            tauSqhat.vsh <- tauSqhat
            tauSqhat.feb <- tauSqhat
            tauSqhat.fes <- tauSqhat
            tauSqhat.reb <- tauSqhat
            tauSqhat.Feb <- tauSqhat
            
            ## CI : tauSq unknown
            CI.bon.unknown   <- CI.bonferroni(sSq,df,xi,theta,sigmaSq,muhat,alpha)
            CI.sSq.unknown   <- CI.t(sSq,df,xi,theta,sigmaSq,muhat,alpha)
            CI.ljs.unknown   <- CI.est(ljs,tauSqhat.ljs,xi,theta,sigmaSq,muhat,alpha)
            CI.opt.unknown   <- CI.est(opt,tauSqhat.opt,xi,theta,sigmaSq,muhat,alpha)
            CI.smy.unknown   <- CI.est(smy,tauSqhat.smy,xi,theta,sigmaSq,muhat,alpha)
            CI.vsh.unknown   <- CI.est(vsh,tauSqhat.vsh,xi,theta,sigmaSq,muhat,alpha)
            CI.feb.unknown   <- CI.est(feb,tauSqhat.feb,xi,theta,sigmaSq,muhat,alpha)
            CI.fes.unknown   <- CI.est(fes,tauSqhat.fes,xi,theta,sigmaSq,muhat,alpha)
            CI.reb.unknown   <- CI.est(reb,tauSqhat.reb,xi,theta,sigmaSq,muhat,alpha)
            CI.Feb.unknown   <- CI.est(Feb,tauSqhat.feb,xi,theta,sigmaSq,muhat,alpha)
            
            ## CI : tauSq unknown - selected
            CI.bon.sel.unknown[i, t,]   <- CI.bon.unknown[max.ind,]
            CI.sSq.sel.unknown[i, t,]   <- CI.sSq.unknown[max.ind,]
            CI.ljs.sel.unknown[i, t,]   <- CI.ljs.unknown[max.ind,]
            CI.opt.sel.unknown[i, t,]   <- CI.opt.unknown[max.ind,]
            CI.smy.sel.unknown[i, t,]   <- CI.smy.unknown[max.ind,]
            CI.vsh.sel.unknown[i, t,]   <- CI.vsh.unknown[max.ind,]
            CI.feb.sel.unknown[i, t,]   <- CI.feb.unknown[max.ind,]
            CI.fes.sel.unknown[i, t,]   <- CI.fes.unknown[max.ind,]
            CI.reb.sel.unknown[i, t,]   <- CI.reb.unknown[max.ind,]
            CI.Feb.sel.unknown[i, t,]   <- CI.Feb.unknown[max.ind,]
          }
      }
        
##        CI.bon.sel.unknown   <-changenames(CI.bon.sel.unknown)
##        CI.sSq.sel.unknown   <-changenames(CI.sSq.sel.unknown)
##        CI.ljs.sel.unknown   <-changenames(CI.ljs.sel.unknown)
##        CI.opt.sel.unknown   <-changenames(CI.opt.sel.unknown)
##        CI.smy.sel.unknown   <-changenames(CI.smy.sel.unknown)
##        CI.vsh.sel.unknown   <-changenames(CI.vsh.sel.unknown)
##        CI.feb.sel.unknown   <-changenames(CI.feb.sel.unknown)
##        CI.fes.sel.unknown   <-changenames(CI.fes.sel.unknown)
##        CI.reb.sel.unknown   <-changenames(CI.reb.sel.unknown)
##        CI.Feb.sel.unknown   <-changenames(CI.Feb.sel.unknown)
        
      
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
        
        
        Result[t,]=c(M[t],bon=bon.sel.unknown.Result,sSq=sSq.sel.unknown.Result,ljs=ljs.sel.unknown.Result,opt=opt.sel.unknown.Result,smy=smy.sel.unknown.Result,vsh=vsh.sel.unknown.Result,feb=feb.sel.unknown.Result,fes=fes.sel.unknown.Result,reb=reb.sel.unknown.Result, Feb=Feb.sel.unknown.Result)
      }
    
      
      ##-- Saving RData -------------------------------------------------------------------------------------------------------------------
      
      path      <- paste("./Data/")
      path1     <- paste("./Figure/")
      settname  <- paste("Indep_prior_",prior,"_df_",df,"_a_",a,"_b_",b,sep="")
      filename  <- paste(path, settname, ".Rdata", sep="")
      save(Result, file = filename)
    }
  
  if(PLOT==TRUE)
    {
      path      <- paste("./Data/")
      path1     <- paste("./Figure/")
      settname  <- paste("Indep_prior_",prior,"_df_",df,"_a_",a,"_b_",b,sep="")
      filename  <- paste(path, settname, ".Rdata", sep="")
      load(filename)
    }

  
  ##-- Saving Figures ----------------------------------------------------------------------------------------------------------------------------
  
  if (prior==1){
    PRIOR<-"Inverse Gamma"
  } else if (prior==2){
    PRIOR<-"Lognormal"
  } else if (prior==3){
    PRIOR<-"Mix INGamma"
  } else if (prior==4){
    PRIOR<-"Mix Lognormal"
  }
  

  par(mar=c(4, 4, 3, 8))
  xloc=-0.2

  path1 = "./Figure/"
  
  XLAB = Result[1:30,1]/(1-Result[1:30,1])
  
  ## Coverage 
  figurename1 <- paste(path1,"Coverage_",settname,".pdf",sep="")
  pdf(figurename1,width=8,height=5)
  par(mar=c(5,5,1,10) )
  plot( XLAB, Result[1:30,2], xlab=TeX(r'($\tau^2$ )'),ylab=TeX(r'(Coverage Prob. of $\theta_{(1)}$)'), ylim=c(0,1), xlim=c(0,max(XLAB)),col="white",main=bquote(paste(.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df) )), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( XLAB,rep(0.95-qnorm(0.975)*sqrt(0.05*0.95/500),30), type="l", col="black",lwd=3,lty=2)
  lines( XLAB,Result[1:30,2] , pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( XLAB,Result[1:30,10], pch=15, type="b", col="blue"  ) # ELJS
  lines( XLAB,Result[1:30,14], pch=7 , type="b", col="purple") # TW
  lines( XLAB,Result[1:30,18], pch=17, type="b", col="green4") # Smyth
  lines( XLAB,Result[1:30,22], pch=16, type="b", col="gold4" ) # Vash
  lines( XLAB,Result[1:30,26], pch=9 , type="b", col="brown" ) # fEBV
  lines( XLAB,Result[1:30,30], pch=8 , type="b", col="orange") # fEBVS
  lines( XLAB,Result[1:30,34], pch=10, type="b", col="green ") # ReBayes
  lines( XLAB,Result[1:30,38], pch=19, type="b", col="red"   ) # FEBV
  par(xpd=TRUE)
  legend("topright",horiz=F, inset = c(-0.3,0.2), legend=c("Nominal level","Bonferroni","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","REBayes","F-EBV"),col=c("black","cyan3","blue","purple","green4","gold4", "brown","orange","green","red"),pch=c(NA,3,15,7,17,16,9,8,10,19), ncol=1, lty=1:1, cex=1)
  dev.off()
  
  ## Length Ratio
  figurename2 <- paste(path1,"LengRatio_",settname,".pdf",sep="")
  pdf(figurename2,width=8,height=5)
  par(mar=c(5,5,1,10) )
  plot( XLAB,Result[1:30,3], xlab=TeX(r'($\tau^2$ )'),ylab="Length Ratio", ylim=c(0,1), xlim=c(0,max(XLAB)),col="white",main=bquote(paste(.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df) )), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( XLAB,Result[1:30,11]/Result[1:30,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( XLAB,Result[1:30,15]/Result[1:30,3], pch=7 , type="b", col="purple") # TW
  lines( XLAB,Result[1:30,19]/Result[1:30,3], pch=17, type="b", col="green4") # Smyth
  lines( XLAB,Result[1:30,23]/Result[1:30,3], pch=16, type="b", col="gold4" ) # Vash
  lines( XLAB,Result[1:30,27]/Result[1:30,3], pch=9 , type="b", col="brown" ) # fEBV
  lines( XLAB,Result[1:30,31]/Result[1:30,3], pch=8 , type="b", col="orange") # fEBVS
  lines( XLAB,Result[1:30,35]/Result[1:30,3], pch=10, type="b", col="green ") # ReBayes
  lines( XLAB,Result[1:30,39]/Result[1:30,3], pch=19, type="b", col="red"   ) # FEBV
  ##  lines( XLAB,Result[1:30,39]/Result[1:30,3], pch=10, type="b", col="green ") # SOC
  ## legend("topright",horiz=F, inset = c(-0.3,0), legend=c("ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","F-EBV","SOC"),col=c("blue","purple","green4","gold4", "brown","orange","red","green"),pch=c(15,7,17,16,9,8,19,10), ncol=1, lty=1:1, cex=1)
  dev.off()


  #-- Showing Figures ----------------------------------------------------------------------------------------------------------------------------
  
  # Coverage 
  
  plot( XLAB,Result[1:30,2],  xlab=TeX(r'($\tau^2$ )'),ylab=TeX(r'(Coverage Prob. of $\theta_{(1)}$)'), ylim=c(0,1), xlim=c(0,max(XLAB)),col="white",main=bquote(paste(.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df) )), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( XLAB,rep(0.95-qnorm(0.975)*sqrt(0.05*0.95/1000),30), type="l", col="black",lwd=3,lty=2)
  lines( XLAB,Result[1:30,2] , pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( XLAB,Result[1:30,10], pch=15, type="b", col="blue"  ) # ELJS
  lines( XLAB,Result[1:30,14], pch=7 , type="b", col="purple") # TW
  lines( XLAB,Result[1:30,18], pch=17, type="b", col="green4") # Smyth
  lines( XLAB,Result[1:30,22], pch=16, type="b", col="gold4" ) # Vash
  lines( XLAB,Result[1:30,26], pch=9 , type="b", col="brown" ) # fEBV
  lines( XLAB,Result[1:30,30], pch=8 , type="b", col="orange") # fEBVS
  lines( XLAB,Result[1:30,34], pch=10, type="b", col="green" ) # Rebayes
  lines( XLAB,Result[1:30,38], pch=19, type="b", col="red"   ) # FEBV
  ## lines( XLAB,Result[1:30,38], pch=10, type="b", col="green ") # SOC
  legend("topright",horiz=F, inset = c(xloc,0), legend=c("Error Bar","Bonferroni","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","REBayes","F-EBV"),col=c("black","cyan3","blue","purple","green4","gold4", "brown","orange","green","red"),pch=c(NA,3,15,7,17,16,9,8,10,19), ncol=1, lty=1:1, cex=1)
  
  # Length Ratio
  plot( XLAB,Result[1:30,3],  xlab=TeX(r'($\tau^2$ )'), ylab="Length Ratio", ylim=c(0,1), xlim=c(0,max(XLAB)),col="white",main=bquote(paste(.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df) )), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( XLAB,Result[1:30,11]/Result[1:30,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( XLAB,Result[1:30,15]/Result[1:30,3], pch=7 , type="b", col="purple") # TW
  lines( XLAB,Result[1:30,19]/Result[1:30,3], pch=17, type="b", col="green4") # Smyth
  lines( XLAB,Result[1:30,23]/Result[1:30,3], pch=16, type="b", col="gold4" ) # Vash
  lines( XLAB,Result[1:30,27]/Result[1:30,3], pch=9 , type="b", col="brown" ) # fEBV
  lines( XLAB,Result[1:30,31]/Result[1:30,3], pch=8 , type="b", col="orange") # fEBVS
  lines( XLAB,Result[1:30,35]/Result[1:30,3], pch=10, type="b", col="green" ) # ReBayes
  lines( XLAB,Result[1:30,39]/Result[1:30,3], pch=19, type="b", col="red"   ) # FEBV
  ## lines( XLAB,Result[1:30,39]/Result[1:30,3], pch=10, type="b", col="green ") # SOC
  legend("topright",horiz=F, inset = c(xloc,0), legend=c("ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","REBayes","F-EBV"),col=c("blue","purple","green4","gold4", "brown","orange","green","red"),pch=c(15,7,17,16,9,8,10,19), ncol=1, lty=1:1, cex=1)

  Result
  
}
