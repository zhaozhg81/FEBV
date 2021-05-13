library( REBayes)
library( Rmosek )

F_variance_Comparison <- function(G, numSim, n, prior, a, b, bw="ucv", PLOT=FALSE)
{
  
  if( PLOT==FALSE){
    
    ## ##-- Setting Parameters ----------------------------------------------------------------------------------------------------------------------------------------------
    
    ## G       <- 1000     # Dimension
    ## numSim  <- 500     # The number of simulation
    ## n       <- 4        # Sample size  (df=n-1)
  
    ##-- Indep Core ------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
    df      <- n-1   # Degree of Freedom
    mu      <- 0     # Mean of prior distrubution
    alpha   <- 0.05
    
    percvec <- c(1/G, 0.01, 0.05)
    N=length(percvec)
    
    GAM     <- 0
    INVGAM  <- 0
    LNORM   <- 0
    MIXIG   <- 0
    MIXGAM  <- 0
    
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
    
    tsse1.sSq <- array(0, c(numSim,N)) 
    tsse1.ljs <- array(0, c(numSim,N))
    tsse1.smy <- array(0, c(numSim,N))
    tsse1.vsh <- array(0, c(numSim,N))
    tsse1.opt <- array(0, c(numSim,N))
    tsse1.feb <- array(0, c(numSim,N))
    tsse1.fes <- array(0, c(numSim,N))
    tsse1.Feb <- array(0, c(numSim,N))
    tsse1.rebayes <- array(0, c(numSim,N))
    
    tsse2.sSq <- array(0, c(numSim,N))
    tsse2.ljs <- array(0, c(numSim,N))
    tsse2.smy <- array(0, c(numSim,N))
    tsse2.vsh <- array(0, c(numSim,N))
    tsse2.opt <- array(0, c(numSim,N))
    tsse2.feb <- array(0, c(numSim,N))
    tsse2.fes <- array(0, c(numSim,N))
    tsse2.Feb <- array(0, c(numSim,N))
    tsse2.rebayes <- array(0, c(numSim,N))
    
  
    tsse0.sSq <- array(0, c(numSim,N))
    tsse0.ljs <- array(0, c(numSim,N))
    tsse0.smy <- array(0, c(numSim,N))
    tsse0.vsh <- array(0, c(numSim,N))
    tsse0.opt <- array(0, c(numSim,N))
    tsse0.feb <- array(0, c(numSim,N))
    tsse0.fes <- array(0, c(numSim,N))
    tsse0.Feb <- array(0, c(numSim,N))
    tsse0.rebayes <- array(0, c(numSim,N))
    
    tssep.sSq <- array(0, c(numSim,N))
    tssep.ljs <- array(0, c(numSim,N))
    tssep.smy <- array(0, c(numSim,N))
    tssep.vsh <- array(0, c(numSim,N))
    tssep.opt <- array(0, c(numSim,N))
    tssep.feb <- array(0, c(numSim,N))
    tssep.fes <- array(0, c(numSim,N))
    tssep.Feb <- array(0, c(numSim,N))
    tssep.rebayes <- array(0, c(numSim,N))
    
    
    
    
    ##-- Simulation ---------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
##    for(t in  1:N){
##      print( paste("t=",t,sep=""))
      for(i in  1:numSim){
          
        if( i%%10==0 ){
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
            bandwidth <- c(0.6,0.6,0.6)
          }
        bon        <- sSq
        sSq        <- sSq
        ljs        <- LJS(sSq,df)
        opt        <- OPT(sSq,df)
        smy        <- SMY(sSq,df)
        vsh        <- VSH(sSq,df)
        feb        <- fEB(sSq,df,bandwidth)
        fes        <- fES(sSq,df,bandwidth)
        Feb        <- FEB(sSq,df)
        rebayes <- GVmix(sSq, array(df,G) )$dy
           
        for(t in 1:N)
        {
          max.ind    <- order(abs(sSq),decreasing=F)[(1):(percvec[t]*(G/1))]
          sel.sSq        <- sSq[max.ind]
          sel.ljs        <- ljs[max.ind]
          sel.opt        <- opt[max.ind]
          sel.smy        <- smy[max.ind]
          sel.vsh        <- vsh[max.ind]
          sel.feb        <- feb[max.ind]
          sel.fes        <- fes[max.ind]
          sel.Feb        <- Feb[max.ind]
          sel.rebayes <- rebayes[max.ind]
          sel.sigmaSq    <- sigmaSq[max.ind]
        
        

        ##--------------------------------------------------------------------------------------------------------------------------------
        ## L0
        tsse0.sSq[i,t] <-(mean((sel.sSq-sel.sigmaSq)^2))  
        tsse0.ljs[i,t] <-(mean((sel.ljs-sel.sigmaSq)^2))  
        tsse0.opt[i,t] <-(mean((sel.opt-sel.sigmaSq)^2))  
        tsse0.smy[i,t] <-(mean((sel.smy-sel.sigmaSq)^2))  
        tsse0.vsh[i,t] <-(mean((sel.vsh-sel.sigmaSq)^2))  
        tsse0.feb[i,t] <-(mean((sel.feb-sel.sigmaSq)^2))  
        tsse0.fes[i,t] <-(mean((sel.fes-sel.sigmaSq)^2))  
        tsse0.Feb[i,t] <-(mean((sel.Feb-sel.sigmaSq)^2))
        tsse0.rebayes[i,t] <-(mean((sel.rebayes-sel.sigmaSq)^2))  
        
        ##L1p
        tssep.sSq[i,t] <-(mean((sel.sSq/sel.sigmaSq-1)^2))
        tssep.ljs[i,t] <-(mean((sel.ljs/sel.sigmaSq-1)^2))
        tssep.opt[i,t] <-(mean((sel.opt/sel.sigmaSq-1)^2))
        tssep.smy[i,t] <-(mean((sel.smy/sel.sigmaSq-1)^2))
        tssep.vsh[i,t] <-(mean((sel.vsh/sel.sigmaSq-1)^2))
        tssep.feb[i,t] <-(mean((sel.feb/sel.sigmaSq-1)^2))
        tssep.fes[i,t] <-(mean((sel.fes/sel.sigmaSq-1)^2))
        tssep.Feb[i,t] <-(mean((sel.Feb/sel.sigmaSq-1)^2))
        tssep.rebayes[i,t] <-(mean((sel.rebayes/sel.sigmaSq-1)^2))
        
        
        ##L1
        tsse1.sSq[i,t] <-(mean((sel.sigmaSq/sel.sSq-1)^2))
        tsse1.ljs[i,t] <-(mean((sel.sigmaSq/sel.ljs-1)^2))
        tsse1.opt[i,t] <-(mean((sel.sigmaSq/sel.opt-1)^2))
        tsse1.smy[i,t] <-(mean((sel.sigmaSq/sel.smy-1)^2))
        tsse1.vsh[i,t] <-(mean((sel.sigmaSq/sel.vsh-1)^2))
        tsse1.feb[i,t] <-(mean((sel.sigmaSq/sel.feb-1)^2))
        tsse1.fes[i,t] <-(mean((sel.sigmaSq/sel.fes-1)^2))
        tsse1.Feb[i,t] <-(mean((sel.sigmaSq/sel.Feb-1)^2))
        tsse1.rebayes[i,t] <-(mean((sel.sigmaSq/sel.rebayes-1)^2))
        
                                        #L2
        tsse2.sSq[i,t] <-(mean((sel.sSq/sel.sigmaSq-log(sel.sSq/sel.sigmaSq)-1)))
        tsse2.ljs[i,t] <-(mean((sel.ljs/sel.sigmaSq-log(sel.ljs/sel.sigmaSq)-1)))
        tsse2.opt[i,t] <-(mean((sel.opt/sel.sigmaSq-log(sel.opt/sel.sigmaSq)-1)))
        tsse2.smy[i,t] <-(mean((sel.smy/sel.sigmaSq-log(sel.smy/sel.sigmaSq)-1)))
        tsse2.vsh[i,t] <-(mean((sel.vsh/sel.sigmaSq-log(sel.vsh/sel.sigmaSq)-1)))
        tsse2.feb[i,t] <-(mean((sel.feb/sel.sigmaSq-log(sel.feb/sel.sigmaSq)-1)))
        tsse2.fes[i,t] <-(mean((sel.fes/sel.sigmaSq-log(sel.fes/sel.sigmaSq)-1)))
        tsse2.Feb[i,t] <-(mean((sel.Feb/sel.sigmaSq-log(sel.Feb/sel.sigmaSq)-1)))
        tsse2.rebayes[i,t] <-(mean((sel.rebayes/sel.sigmaSq-log(sel.rebayes/sel.sigmaSq)-1)))
        
      }
      
    }
    
    tmse0.sSq<-apply(tsse0.sSq,2,mean);  tmse1.sSq<-apply(tsse1.sSq,2,mean); tmse2.sSq<-apply(tsse2.sSq,2,mean);  tmsep.sSq<-apply(tssep.sSq,2,mean);
    tmse0.ljs<-apply(tsse0.ljs,2,mean);  tmse1.ljs<-apply(tsse1.ljs,2,mean); tmse2.ljs<-apply(tsse2.ljs,2,mean);  tmsep.ljs<-apply(tssep.ljs,2,mean); 
    tmse0.opt<-apply(tsse0.opt,2,mean);  tmse1.opt<-apply(tsse1.opt,2,mean); tmse2.opt<-apply(tsse2.opt,2,mean);  tmsep.opt<-apply(tssep.opt,2,mean); 
    tmse0.smy<-apply(tsse0.smy,2,mean);  tmse1.smy<-apply(tsse1.smy,2,mean); tmse2.smy<-apply(tsse2.smy,2,mean);  tmsep.smy<-apply(tssep.smy,2,mean); 
    tmse0.vsh<-apply(tsse0.vsh,2,mean);  tmse1.vsh<-apply(tsse1.vsh,2,mean); tmse2.vsh<-apply(tsse2.vsh,2,mean);  tmsep.vsh<-apply(tssep.vsh,2,mean); 
    tmse0.feb<-apply(tsse0.feb,2,mean);  tmse1.feb<-apply(tsse1.feb,2,mean); tmse2.feb<-apply(tsse2.feb,2,mean);  tmsep.feb<-apply(tssep.feb,2,mean); 
    tmse0.fes<-apply(tsse0.fes,2,mean);  tmse1.fes<-apply(tsse1.fes,2,mean); tmse2.fes<-apply(tsse2.fes,2,mean);  tmsep.fes<-apply(tssep.fes,2,mean); 
    tmse0.Feb<-apply(tsse0.Feb,2,mean);  tmse1.Feb<-apply(tsse1.Feb,2,mean); tmse2.Feb<-apply(tsse2.Feb,2,mean);  tmsep.Feb<-apply(tssep.Feb,2,mean);
    tmse0.rebayes<-apply(tsse0.rebayes,2,mean);  tmse1.rebayes<-apply(tsse1.rebayes,2,mean); tmse2.rebayes<-apply(tsse2.rebayes,2,mean);  tmsep.rebayes<-apply(tssep.rebayes,2,mean); 
    
    
    tmsd0.sSq<-apply(tsse0.sSq,2,sd);  tmsd1.sSq<-apply(tsse1.sSq,2,sd); tmsd2.sSq<-apply(tsse2.sSq,2,sd);  tmsdp.sSq<-apply(tssep.sSq,2,sd);
    tmsd0.ljs<-apply(tsse0.ljs,2,sd);  tmsd1.ljs<-apply(tsse1.ljs,2,sd); tmsd2.ljs<-apply(tsse2.ljs,2,sd);  tmsdp.ljs<-apply(tssep.ljs,2,sd); 
    tmsd0.opt<-apply(tsse0.opt,2,sd);  tmsd1.opt<-apply(tsse1.opt,2,sd); tmsd2.opt<-apply(tsse2.opt,2,sd);  tmsdp.opt<-apply(tssep.opt,2,sd); 
    tmsd0.smy<-apply(tsse0.smy,2,sd);  tmsd1.smy<-apply(tsse1.smy,2,sd); tmsd2.smy<-apply(tsse2.smy,2,sd);  tmsdp.smy<-apply(tssep.smy,2,sd); 
    tmsd0.vsh<-apply(tsse0.vsh,2,sd);  tmsd1.vsh<-apply(tsse1.vsh,2,sd); tmsd2.vsh<-apply(tsse2.vsh,2,sd);  tmsdp.vsh<-apply(tssep.vsh,2,sd); 
    tmsd0.feb<-apply(tsse0.feb,2,sd);  tmsd1.feb<-apply(tsse1.feb,2,sd); tmsd2.feb<-apply(tsse2.feb,2,sd);  tmsdp.feb<-apply(tssep.feb,2,sd); 
    tmsd0.fes<-apply(tsse0.fes,2,sd);  tmsd1.fes<-apply(tsse1.fes,2,sd); tmsd2.fes<-apply(tsse2.fes,2,sd);  tmsdp.fes<-apply(tssep.fes,2,sd); 
    tmsd0.Feb<-apply(tsse0.Feb,2,sd);  tmsd1.Feb<-apply(tsse1.Feb,2,sd); tmsd2.Feb<-apply(tsse2.Feb,2,sd);  tmsdp.Feb<-apply(tssep.Feb,2,sd); 
    tmsd0.rebayes<-apply(tsse0.rebayes,2,sd);  tmsd1.rebayes<-apply(tsse1.rebayes,2,sd); tmsd2.rebayes<-apply(tsse2.rebayes,2,sd);  tmsdp.rebayes<-apply(tssep.rebayes,2,sd); 
    
    
    
    L0<-round(cbind(tmse0.sSq, tmse0.ljs, tmse0.opt, tmse0.smy, tmse0.vsh, tmse0.feb, tmse0.fes, tmse0.rebayes, tmse0.Feb ),3)
    L1<-round(cbind(tmse1.sSq, tmse1.ljs, tmse1.opt, tmse1.smy, tmse1.vsh, tmse1.feb, tmse1.fes, tmse1.rebayes, tmse1.Feb),3)
    L2<-round(cbind(tmse2.sSq, tmse2.ljs, tmse2.opt, tmse2.smy, tmse2.vsh, tmse2.feb, tmse2.fes, tmse2.rebayes, tmse2.Feb),3)
    Lp<-round(cbind(tmsep.sSq, tmsep.ljs, tmsep.opt, tmsep.smy, tmsep.vsh, tmsep.feb, tmsep.fes, tmsep.rebayes, tmsep.Feb),3)
  
    Lsd0<-round(cbind(tmsd0.sSq, tmsd0.ljs, tmsd0.opt, tmsd0.smy, tmsd0.vsh, tmsd0.feb,  tmsd0.fes, tmsd0.rebayes, tmsd0.Feb),3)
    Lsd1<-round(cbind(tmsd1.sSq, tmsd1.ljs, tmsd1.opt, tmsd1.smy, tmsd1.vsh, tmsd1.feb,  tmsd1.fes, tmsd1.rebayes, tmsd1.Feb),3)
    Lsd2<-round(cbind(tmsd2.sSq, tmsd2.ljs, tmsd2.opt, tmsd2.smy, tmsd2.vsh, tmsd2.feb,  tmsd2.fes, tmsd2.rebayes, tmsd2.Feb),3)
    Lsdp<-round(cbind(tmsdp.sSq, tmsdp.ljs, tmsdp.opt, tmsdp.smy, tmsdp.vsh, tmsdp.feb,  tmsdp.fes, tmsdp.rebayes, tmsdp.Feb),3)
    
    ##-- Saving RData -------------------------------------------------------------------------------------------------------------------
    
    path      <- paste("./Data/")
    path1     <- paste("./Figure/")
    
    settname0  <- paste("VarComp_L0_prior_",prior,"_G_",G,"_a_",a,"_b_",b,sep="")
    filename   <- paste(path, settname0, ".Rdata", sep="")
    save(L0, file = filename)
    
    settname1  <- paste("VarComp_L1_prior_",prior,"_G_",G,"_a_",a,"_b_",b,sep="")
    filename   <- paste(path, settname1, ".Rdata", sep="")
    save(L1, file = filename)
    
    settname2  <- paste("VarComp_L2_prior_",prior,"_G_",G,"_a_",a,"_b_",b,sep="")
    filename   <- paste(path, settname2, ".Rdata", sep="")
    save(L2, file = filename)
    
    settnamep  <- paste("VarComp_Lp_prior_",prior,"_G_",G,"_a_",a,"_b_",b,sep="")
    filename   <- paste(path, settnamep, ".Rdata", sep="")
    save(Lp, file = filename)
  }
  
  if( PLOT==TRUE)
    {
      path      <- paste("./Data/")
      path1     <- paste("./Figure/")
      
      settname0  <- paste("VarComp_L0_prior_",prior,"_G_",G,"_a_",a,"_b_",b,sep="")
      filename   <- paste(path, settname0, ".Rdata", sep="")
      load(filename)
      
      settname1  <- paste("VarComp_L1_prior_",prior,"_G_",G,"_a_",a,"_b_",b,sep="")
      filename   <- paste(path, settname1, ".Rdata", sep="")
      load(filename)
      
      settname2  <- paste("VarComp_L2_prior_",prior,"_G_",G,"_a_",a,"_b_",b,sep="")
      filename   <- paste(path, settname2, ".Rdata", sep="")
      load(filename)
      
      settnamep  <- paste("VarComp_Lp_prior_",prior,"_G_",G,"_a_",a,"_b_",b,sep="")
      filename   <- paste(path, settnamep, ".Rdata", sep="")
      load(filename)
    }

  ##-- Saving Figures -----------------------------------------------------------------------------------------------------------------
  
  path1     <- paste("./Figure/")
  
  ##-- Saving Figures ----------------------------------------------------------------------------------------------------------------------------
  
  if (prior==1){
    PRIOR<-"InGamma"
  } else if (prior==2){
    PRIOR<-"Lognormal"
  } else if (prior==4){
    PRIOR<-"Mix Lognormal"
  } else if (prior==3){
    PRIOR<-"Mix InGamma"
  }
  
  par(mfrow=c(2,2))
  namelist=c("sSq",	"ELJS",	"TW",	"Smyth",	"Vash",	"fEBV",	"fEBVS", "REBayes",	"FEBV")
  
  figurename1 <- paste(path1,"Risk_var_",settname0,"_max_five_percent.pdf",sep="")
  pdf(figurename1,width=12,height=10)
  par(mfrow=c(1,2))
  barplot(log(L0[1,],10), ylab="log(MSE)", main=bquote(paste("[L0, Max]  ",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
  barplot(log(L0[3,],10), ylab="log(MSE)", main=bquote(paste("[L0, 5%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
  dev.off()
  
  figurename1 <- paste(path1,"Risk_var_",settname1,"_max_five_percent.pdf",sep="")
  pdf(figurename1,width=12,height=10)
  par(mfrow=c(1,2))
  barplot(log(L1[1,],10), ylab="log(MSE)", main=bquote(paste("[L1, Max]  ",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
  barplot(log(L1[3,],10), ylab="log(MSE)", main=bquote(paste("[L1, 5%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
  dev.off()
  
  figurename1 <- paste(path1,"Risk_var_",settname2,"_max_five_percent.pdf",sep="")
  pdf(figurename1,width=12,height=10)
  par(mfrow=c(1,2))
  barplot(log(L2[1,],10), ylab="log(MSE)", main=bquote(paste("[L2, Max]  ",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
  barplot(log(L2[3,],10), ylab="log(MSE)", main=bquote(paste("[L2, 5%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
  dev.off()
  
  figurename1 <- paste(path1,"Risk_var_",settnamep,"_max_five_percent.pdf",sep="")
  pdf(figurename1,width=12,height=10)
  par(mfrow=c(1,2))
  barplot(log(Lp[1,],10), ylab="log(MSE)", main=bquote(paste("[L1p, Max]  ",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
  barplot(log(Lp[3,],10), ylab="log(MSE)", main=bquote(paste("[L1p, 5%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
  dev.off()

  

}
