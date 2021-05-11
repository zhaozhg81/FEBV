library(latex2exp)
library(REBayes)


RANK <- function( xi, sigmaSqHat )
{
  p = length(xi)
  T.stat = xi/sqrt( sigmaSqHat)
  
  rank = order(abs(T.stat), decreasing=TRUE)
  rank
}


F_Testing <- function(G, numSim, N, n, prior, a, b, perc=0.001, bw='ucv', PLOT=FALSE, REJ=100)
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
    
    TP.sSq = array(0, c(numSim, N) )
    TP.ljs = TP.sSq 
    TP.opt = TP.sSq 
    TP.smy = TP.sSq 
    TP.vsh = TP.sSq 
    TP.feb = TP.sSq 
    TP.fes = TP.sSq 
    TP.rebayes = TP.sSq 
    TP.Feb = TP.sSq 
    
    

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
      ##vsh        <- VSH(sSq,df)
      opt        <- OPT(sSq,df)
      ##feb        <- fEB(sSq,df,bandwidth)
      ##fes        <- fES(sSq,df,bandwidth)
      Feb        <- FEB(sSq,df)
      rebayes    <- GVmix(sSq, array(df,G) )$dy
      vsh = Feb
      feb = Feb
      fes = Feb
      
      for(t in 1:N)
      {
        tauSq      <- M[t]/(1-M[t])
        tmp <- runif(G)
        
        true.ind <- (tmp<0.1)
        theta      <- (tmp<0.1)*rnorm(G,mu,sqrt(tauSq))
        
        xi         <- rnorm(G,theta,sqrt(sigmaSq))
        
        rank.sSq <- RANK( xi, sSq )
        rank.ljs <- RANK( xi, ljs )
        rank.opt <- RANK( xi, opt )
        rank.smy <- RANK( xi, smy )
        rank.vsh <- RANK( xi, vsh )
        rank.feb <- RANK( xi, feb )
        rank.fes <- RANK( xi, fes )
        rank.rebayes <- RANK( xi, rebayes )
        rank.Feb <- RANK( xi, Feb )
        
        TP.sSq[i, t] =  sum( true.ind[ rank.sSq[1:REJ] ] )
        TP.ljs[i, t] =  sum( true.ind[ rank.ljs[1:REJ] ] )
        TP.opt[i, t] =  sum( true.ind[ rank.opt[1:REJ] ] )
        TP.smy[i, t] =  sum( true.ind[ rank.smy[1:REJ] ] )
        TP.vsh[i, t] =  sum( true.ind[ rank.vsh[1:REJ] ] )
        TP.feb[i, t] =  sum( true.ind[ rank.feb[1:REJ] ] )
        TP.fes[i, t] =  sum( true.ind[ rank.fes[1:REJ] ] )
        TP.rebayes[i, t] =  sum( true.ind[ rank.rebayes[1:REJ] ] )
        TP.Feb[i, t] =  sum( true.ind[ rank.Feb[1:REJ] ] )
        
        
      }
    }
    

    sSq.res = apply( TP.sSq, 2, mean)
    ljs.res = apply( TP.ljs, 2, mean)
    opt.res = apply( TP.opt, 2, mean)
    smy.res = apply( TP.smy, 2, mean)
    vsh.res = apply( TP.vsh, 2, mean)
    feb.res = apply( TP.feb, 2, mean)
    fes.res = apply( TP.fes, 2, mean)
    rebayes.res = apply( TP.rebayes, 2, mean)
    Feb.res = apply( TP.Feb, 2, mean)
    
    Result <-cbind(M, sSq.res, ljs.res, opt.res, smy.res, vsh.res, feb.res, fes.res, rebayes.res, Feb.res)
    
    ##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ##-- Saving RData -------------------------------------------------------------------------------------------------------------------
    
    path      <- paste("./Data/")
    settname  <- paste("Testing_prior_",prior,"_G_",G,"_df_",df,"_a_",a,"_b_",b,"_perc_",perc,sep="")
    filename  <- paste(path, settname, ".Rdata", sep="")
    save(Result, file = filename)
    
  }
  
  if(PLOT==TRUE)
  {
    ##-- Saving RData -------------------------------------------------------------------------------------------------------------------
    
    path      <- paste("./Data/")
    settname  <- paste("Testing_prior_",prior,"_G_",G,"_df_",df,"_a_",a,"_b_",b,"_perc_",perc,sep="")
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
  
  ## Comparison of Estimators for theta
  figurename1 <- paste(path1,"Testing_",settname,".pdf",sep="")
  pdf(figurename1,width=8,height=5)
  par(mar=c(5,6,1,10) )
  plot( Result[,1],Result[,2],xlab=TeX(r'($\tau^2/(\tau^2+E\sigma^2)$ )'),ylab=TeX(r'(MSE of $\theta_{(1)}$ )'), ylim=c(0,max(Result[,2])), xlim=c(0,1),col="white",main=bquote(paste(.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df) )), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( Result[,1],Result[,2], pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( Result[,1],Result[,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result[,1],Result[,4], pch=7 , type="b", col="purple") # TW
  lines( Result[,1],Result[,5], pch=17, type="b", col="green4") # Smyth
  lines( Result[,1],Result[,6], pch=16, type="b", col="gold4" ) # Vash
  lines( Result[,1],Result[,7], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result[,1],Result[,8], pch=8 , type="b", col="orange") # fEBVS
  lines( Result[,1],Result[,9], pch=10, type="b", col="green ") # ReBayes
  lines( Result[,1],Result[,10], pch=19, type="b", col="red"  ) # FEBV
  
  if(prior <=3 )
  {
    par(xpd=TRUE)
    legend("topright",horiz=F, inset = c(-0.3,0.2), legend=c("sSq","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS", "REBayes","F-EBV"),col=c("cyan3","blue","purple","green4","gold4", "brown","orange","green","red"),pch=c(3,15,7,17,16,9,8,10, 19), ncol=1, lty=1:1, cex=1)
  }
  dev.off()
  
  path1     <-"/home/zhaozhg/Dropbox/Apps/Overleaf/On F-modelling based Empirical Bayes Estimation of Variances/figure/"
  
  ## Comparison of Estimators for theta
  figurename1 <- paste(path1,"Testing_",settname,".pdf",sep="")
  pdf(figurename1,width=8,height=5)
  par(mar=c(5,6,1,10) )
  plot( Result[,1],Result[,2],xlab=TeX(r'($\tau^2/(\tau^2+E\sigma^2)$ )'),ylab=TeX(r'(MSE of $\theta_{(1)}$ )'), ylim=c(0,max(Result[,2])), xlim=c(0,1),col="white",main=bquote(paste(.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df) )), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( Result[,1],Result[,2], pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( Result[,1],Result[,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result[,1],Result[,4], pch=7 , type="b", col="purple") # TW
  lines( Result[,1],Result[,5], pch=17, type="b", col="green4") # Smyth
  lines( Result[,1],Result[,6], pch=16, type="b", col="gold4" ) # Vash
  lines( Result[,1],Result[,7], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result[,1],Result[,8], pch=8 , type="b", col="orange") # fEBVS
  lines( Result[,1],Result[,9], pch=10, type="b", col="green ") # ReBayes
  lines( Result[,1],Result[,10], pch=19, type="b", col="red"  ) # FEBV
  
  if(prior <=3 )
  {
    par(xpd=TRUE)
    legend("topright",horiz=F, inset = c(-0.3,0.2), legend=c("sSq","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS", "REBayes","F-EBV"),col=c("cyan3","blue","purple","green4","gold4", "brown","orange","green","red"),pch=c(3,15,7,17,16,9,8,10, 19), ncol=1, lty=1:1, cex=1)
  }
  dev.off()  
  
  
  ##-- Showing Figures ----------------------------------------------------------------------------------------------------------------------------
  par(mar=c(5,6,1,10) )
  plot( Result[,1],Result[,2],xlab=TeX(r'($\tau^2/(\tau^2+E\sigma^2)$ )'),ylab=TeX(r'(MSE of $\theta_{(1)}$ )'), ylim=c(0,max(Result[,2])), xlim=c(0,1),col="white",main=bquote(paste(.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df) )), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( Result[,1],Result[,2], pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( Result[,1],Result[,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result[,1],Result[,4], pch=7 , type="b", col="purple") # TW
  lines( Result[,1],Result[,5], pch=17, type="b", col="green4") # Smyth
  lines( Result[,1],Result[,6], pch=16, type="b", col="gold4" ) # Vash
  lines( Result[,1],Result[,7], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result[,1],Result[,8], pch=8 , type="b", col="orange") # fEBVS
  lines( Result[,1],Result[,9], pch=10, type="b", col="green ") # ReBayes
  lines( Result[,1],Result[,10], pch=19, type="b", col="red"  ) # FEBV
  
  if(prior <=3 )
  {
    par(xpd=TRUE)
    legend("topright",horiz=F, inset = c(-0.3,0.2), legend=c("sSq","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS", "REBayes","F-EBV"),col=c("cyan3","blue","purple","green4","gold4", "brown","orange","green","red"),pch=c(3,15,7,17,16,9,8,10, 19), ncol=1, lty=1:1, cex=1)
  }
  dev.off()
  
  Result
}
