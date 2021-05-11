library(latex2exp)
library(REBayes)

F_mean_Comparison <- function(G, numSim, N, n, prior, a, b, perc=0.001, bw='ucv', PLOT=FALSE)
{
  
  if(PLOT==FALSE)
    {
      
      ## G       <- 1000   # Dimension
      ## numSim  <- 500  # The number of simulation
      ## N       <- 20     # Number of ticks for X-axis
      ## n       <- 4      # Sample size  (df=n-1)

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
      smy     <- array(0, G)
      vsh     <- array(0, G)
      feb     <- array(0, G)
      fes     <- array(0, G)
      Feb     <- array(0, G)
      rebayes <- array(0, G)
      
      sel.sSq.m       <- array(0, c(numSim,N))
      sel.ljs.m       <- sel.sSq.m
      sel.smy.m       <- sel.sSq.m
      sel.vsh.m       <- sel.sSq.m
      sel.opt.m       <- sel.sSq.m
      sel.feb.m       <- sel.sSq.m
      sel.fes.m       <- sel.sSq.m
      sel.Feb.m       <- sel.sSq.m
      sel.rebayes.m   <- sel.sSq.m
      

      Result           <- array(0, c(N,33))
      
      M=seq(0,1,0.999/N)
      M=M[2:(N+1)]
      
      ##-- Simulation ---------------------------------------------------------------------------------------------------------------------------------------------------------
        ## for(t in  1:N){
        ## print( paste("t=",t,sep=""))
      for(i in  1:numSim){
          if( i%% 20==0 ){
            print(date())
            print( paste( "number of simulation:", i, sep="") )
          }

        sigmaSq <- GenerateVar(G, prior, a, b)
          
          sSq        <- sigmaSq*rchisq(G,df)/df
          
          

          if( bw=='ucv')
            {
              if( i==1 )
                {
                  bandwidth <- array(0, 3)
                  bandwidth[1] <- h.ucv(sSq, deriv.order=0)$h
                  bandwidth[2] <- h.ucv(sSq, deriv.order=1)$h
                  bandwidth[3] <- h.ucv(sSq, deriv.order=2)$h                
                }
            }else{
              bandwidth <- array(0.6, 3)
            }
          
          bon        <- sSq
          sSq        <- sSq
          ljs        <- LJS(sSq,df)
          smy        <- SMY(sSq,df)
        vsh        <- VSH(sSq,df)
          opt        <- OPT(sSq,df)
          feb        <- fEB(sSq,df,bandwidth)
          fes        <- fES(sSq,df,bandwidth)
          Feb        <- FEB(sSq,df)
          rebayes    <- GVmix(sSq, array(df,G) )$dy
          
        for(t in 1:N)
        {
          tauSq      <- M[t]/(1-M[t])
          theta      <- rnorm(G,mu,sqrt(tauSq))
          xi         <- rnorm(G,theta,sqrt(sigmaSq))
          
          muhat      <- mean(xi) #sum((xi/sSq))/sum(1/sSq)
          T.stat     <- xi/sqrt(sSq)
          max.ind    <- order(abs(T.stat),decreasing=TRUE)[((G/G):(perc*G/1))]
          
          
          tauSqhat.sSq<- tauSqfn(sSq,xi,muhat,trunc=TRUE)
          tauSqhat.ljs<- tauSqfn(ljs,xi,muhat,trunc=TRUE)
          tauSqhat.smy<- tauSqfn(smy,xi,muhat,trunc=TRUE)
          tauSqhat.vsh<- tauSqfn(vsh,xi,muhat,trunc=TRUE)
          tauSqhat.opt<- tauSqfn(opt,xi,muhat,trunc=TRUE)
          tauSqhat.fes<- tauSqfn(fes,xi,muhat,trunc=TRUE)
          tauSqhat.feb<- tauSqfn(feb,xi,muhat,trunc=TRUE)
          tauSqhat.Feb<- tauSqfn(Feb,xi,muhat,trunc=TRUE)
          tauSqhat.rebayes<- tauSqfn(rebayes,xi,muhat,trunc=TRUE)

          tauSqhat<- tauSqfn(sSq,xi,muhat,trunc=TRUE)
          
          tauSqhat.sSq<- tauSqhat
          tauSqhat.ljs<- tauSqhat
          tauSqhat.smy<- tauSqhat
          tauSqhat.vsh<- tauSqhat
          tauSqhat.opt<- tauSqhat
          tauSqhat.fes<- tauSqhat
          tauSqhat.feb<- tauSqhat
          tauSqhat.Feb<- tauSqhat
          tauSqhat.rebayes<- tauSqhat
          
          ## tauSqhat.sSq[i]<- tauSq
          ## tauSqhat.ljs[i]<- tauSq
          ## tauSqhat.smy[i]<- tauSq
          ## tauSqhat.vsh[i]<- tauSq
          ## tauSqhat.opt[i]<- tauSq
          ## tauSqhat.fes[i]<- tauSq
          ## tauSqhat.feb[i]<- tauSq
          ## tauSqhat.Feb[i]<- tauSq
          ## tauSqhat.rebayes[i]<- tauSq

          
          ## CI : tauSq unknown
          
          sSq.m  <- (tauSqhat.sSq)/(sSq+tauSqhat.sSq)*xi+(sSq)/(sSq+tauSqhat.sSq)*muhat
          ljs.m  <- (tauSqhat.ljs)/(ljs+tauSqhat.ljs)*xi+(ljs)/(ljs+tauSqhat.ljs)*muhat
          opt.m  <- (tauSqhat.opt)/(opt+tauSqhat.opt)*xi+(opt)/(opt+tauSqhat.opt)*muhat
          smy.m  <- (tauSqhat.smy)/(smy+tauSqhat.smy)*xi+(smy)/(smy+tauSqhat.smy)*muhat
          vsh.m  <- (tauSqhat.vsh)/(vsh+tauSqhat.vsh)*xi+(vsh)/(vsh+tauSqhat.vsh)*muhat
          feb.m  <- (tauSqhat.feb)/(feb+tauSqhat.feb)*xi+(feb)/(feb+tauSqhat.feb)*muhat
          fes.m  <- (tauSqhat.fes)/(feb+tauSqhat.fes)*xi+(fes)/(fes+tauSqhat.fes)*muhat
          Feb.m  <- (tauSqhat.Feb)/(Feb+tauSqhat.Feb)*xi+(Feb)/(Feb+tauSqhat.Feb)*muhat
          rebayes.m  <- (tauSqhat.rebayes)/(rebayes+tauSqhat.rebayes)*xi+(rebayes)/(rebayes+tauSqhat.rebayes)*muhat
          
          dif.sSq.m <- (sSq.m-theta)^2
          dif.ljs.m <- (ljs.m-theta)^2
          dif.opt.m <- (opt.m-theta)^2
          dif.smy.m <- (smy.m-theta)^2
          dif.vsh.m <- (vsh.m-theta)^2
          dif.feb.m <- (feb.m-theta)^2
          dif.fes.m <- (feb.m-theta)^2
          dif.Feb.m <- (Feb.m-theta)^2
          dif.rebayes.m <- (rebayes.m-theta)^2
          
          
          ## CI : tauSq unknown - selected
          sel.sSq.m[i,t]   <- mean(dif.sSq.m[max.ind])
          sel.ljs.m[i,t]   <- mean(dif.ljs.m[max.ind])
          sel.opt.m[i,t]   <- mean(dif.opt.m[max.ind])
          sel.smy.m[i,t]   <- mean(dif.smy.m[max.ind])
          sel.vsh.m[i,t]   <- mean(dif.vsh.m[max.ind])
          sel.feb.m[i,t]   <- mean(dif.feb.m[max.ind])
          sel.fes.m[i,t]   <- mean(dif.fes.m[max.ind])
          sel.Feb.m[i,t]   <- mean(dif.Feb.m[max.ind])
          sel.rebayes.m[i,t]   <- mean(dif.rebayes.m[max.ind])
       }
      }

      
        sSqR=apply(sel.sSq.m,2,mean)#/mean(sel.Feb.m);
        ljsR=apply(sel.ljs.m,2,mean)#/mean(sel.Feb.m);
        optR=apply(sel.opt.m,2,mean)#/mean(sel.Feb.m);
        smyR=apply(sel.smy.m,2,mean)#/mean(sel.Feb.m);
        vshR=apply(sel.vsh.m,2,mean)#/mean(sel.Feb.m);
        febR=apply(sel.feb.m,2,mean)#/mean(sel.Feb.m);
        fesR=apply(sel.fes.m,2,mean)#/mean(sel.Feb.m);
        FebR=apply(sel.Feb.m,2,mean)#/mean(sel.Feb.m);
        rebayesR=apply(sel.rebayes.m,2,mean)#/mean(sel.Feb.m);
        
        
 
      Result <-cbind(M, sSqR, ljsR, optR, smyR, vshR, febR, fesR, rebayesR, FebR)
      
      ##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      ##-- Saving RData -------------------------------------------------------------------------------------------------------------------

      path      <- paste("./Data/")
      settname  <- paste("MeanComp_prior_",prior,"_G_",G,"_df_",df,"_a_",a,"_b_",b,"_perc_",perc,sep="")
      filename  <- paste(path, settname, ".Rdata", sep="")
      save(Result, file = filename)
      
    }
  
  if(PLOT==TRUE)
    {
      ##-- Saving RData -------------------------------------------------------------------------------------------------------------------
      
      path      <- paste("./Data/")
      settname  <- paste("MeanComp_prior_",prior,"_G_",G,"_df_",df,"_a_",a,"_b_",b,"_perc_",perc,sep="")
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
  
  par(mar=c(4, 5, 3, 8))
  xloc=-0.2
  path1     <-"./Figure/"
  
  XLAB = Result[1:30,1]/(1-Result[1:30,1])
  
  ## Comparison of Estimators for theta
  figurename1 <- paste(path1,"MSE_mean_comp_",settname,".pdf",sep="")
  pdf(figurename1,width=8,height=5)
  par(mar=c(5,6,1,10) )
  plot( XLAB,Result[1:30,2],xlab=TeX(r'($\tau^2$ )'),ylab=TeX(r'(MSE of $\theta_{(1)}$ )'), ylim=c(0,max(Result[1:30,2])), xlim=c(0,max(XLAB)),col="white",main=bquote(paste(.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df) )), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( XLAB,Result[1:30,2], pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( XLAB,Result[1:30,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( XLAB,Result[1:30,4], pch=7 , type="b", col="purple") # TW
  lines( XLAB,Result[1:30,5], pch=17, type="b", col="green4") # Smyth
  lines( XLAB,Result[1:30,6], pch=16, type="b", col="gold4" ) # Vash
  lines( XLAB,Result[1:30,7], pch=9 , type="b", col="brown" ) # fEBV
  lines( XLAB,Result[1:30,8], pch=8 , type="b", col="orange") # fEBVS
  lines( XLAB,Result[1:30,9], pch=10, type="b", col="green ") # ReBayes
  lines( XLAB,Result[1:30,10], pch=19, type="b", col="red"  ) # FEBV

  if(prior <= 2 )
    {
      par(xpd=TRUE)
      legend("topright",horiz=F, inset = c(-0.3,0.2), legend=c("sSq","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS", "REBayes","F-EBV"),col=c("cyan3","blue","purple","green4","gold4", "brown","orange","green","red"),pch=c(3,15,7,17,16,9,8,10, 19), ncol=1, lty=1:1, cex=1)
    }
  dev.off()

  path1     <-"/home/zhaozhg/Dropbox/Apps/Overleaf/On F-modelling based Empirical Bayes Estimation of Variances/figure/"
  
  ## Comparison of Estimators for theta
  figurename1 <- paste(path1,"MSE_mean_comp_",settname,".pdf",sep="")
  pdf(figurename1,width=8,height=5)
  par(mar=c(5,6,1,10) )
  plot( XLAB,Result[1:30,2],xlab=TeX(r'($\tau^2$ )'),ylab=TeX(r'(MSE of $\theta_{(1)}$ )'), ylim=c(0,max(Result[1:30,2])), xlim=c(0,max(XLAB)),col="white",main=bquote(paste(.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df) )), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( XLAB,Result[1:30,2], pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( XLAB,Result[1:30,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( XLAB,Result[1:30,4], pch=7 , type="b", col="purple") # TW
  lines( XLAB,Result[1:30,5], pch=17, type="b", col="green4") # Smyth
  lines( XLAB,Result[1:30,6], pch=16, type="b", col="gold4" ) # Vash
  lines( XLAB,Result[1:30,7], pch=9 , type="b", col="brown" ) # fEBV
  lines( XLAB,Result[1:30,8], pch=8 , type="b", col="orange") # fEBVS
  lines( XLAB,Result[1:30,9], pch=10, type="b", col="green ") # ReBayes
  lines( XLAB,Result[1:30,10], pch=19, type="b", col="red"  ) # FEBV

  if(prior <= 2)
    {
      par(xpd=TRUE)
      legend("topright",horiz=F, inset = c(-0.3,0.2), legend=c("sSq","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS", "REBayes","F-EBV"),col=c("cyan3","blue","purple","green4","gold4", "brown","orange","green","red"),pch=c(3,15,7,17,16,9,8,10, 19), ncol=1, lty=1:1, cex=1)
    }
  dev.off()  

    
  ##-- Showing Figures ----------------------------------------------------------------------------------------------------------------------------
  par(mar=c(5,6,1,10) )
  plot( XLAB,Result[1:30,2],xlab=TeX(r'($\tau^2$ )'),ylab=TeX(r'(MSE of $\theta_{(1)}$ )'), ylim=c(0,max(Result[1:30,2])), xlim=c(0,max(XLAB)),col="white",main=bquote(paste(.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df) )), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( XLAB,Result[1:30,2], pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( XLAB,Result[1:30,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( XLAB,Result[1:30,4], pch=7 , type="b", col="purple") # TW
  lines( XLAB,Result[1:30,5], pch=17, type="b", col="green4") # Smyth
  lines( XLAB,Result[1:30,6], pch=16, type="b", col="gold4" ) # Vash
  lines( XLAB,Result[1:30,7], pch=9 , type="b", col="brown" ) # fEBV
  lines( XLAB,Result[1:30,8], pch=8 , type="b", col="orange") # fEBVS
  lines( XLAB,Result[1:30,9], pch=10, type="b", col="green ") # ReBayes
  lines( XLAB,Result[1:30,10], pch=19, type="b", col="red"  ) # FEBV

  if(prior <=2 )
    {
      par(xpd=TRUE)
      legend("topright",horiz=F, inset = c(-0.3,0.2), legend=c("sSq","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS", "REBayes","F-EBV"),col=c("cyan3","blue","purple","green4","gold4", "brown","orange","green","red"),pch=c(3,15,7,17,16,9,8,10, 19), ncol=1, lty=1:1, cex=1)
    }
  dev.off()

  Result
}
