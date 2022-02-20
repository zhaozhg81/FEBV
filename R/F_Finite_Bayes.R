library(REBayes)
library(Rmosek)

F_FinBayes<- function(G, numSim, n, prior, a, b, PLOT=FALSE, bw='ucv')
{

  

  if(PLOT==FALSE)
  {
    
    ## ##-- Setting Parameters ----------------------------------------------------------------------------------------------------------------------------------------------
    
    ## G       <- 1000   # Dimension
    ## numSim  <- 500    # The number of simulation
    ## N       <- 1      # Number of ticks for X-axis
    ## n       <- 4      # Sample size  (df=n-1)
    
    ##-- Indep Core ------------------------------------------------------------------------------------------------------------------------------------------------------
    N <- 1
    df      <- n-1   # Degree of Freedom
    mu      <- 0     # Mean of prior distrubution
    alpha   <- 0.05
    
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
    opt     <- array(0, G)
    smy     <- array(0, G)
    modified.smy     <- array(0, G)
    vsh     <- array(0, G)
    feb     <- array(0, G)
    fes     <- array(0, G)
    Feb     <- array(0, G)
    rebayes <- array(0, G)
    
    sSq0     <- array(0, numSim)
    ljs0     <- array(0, numSim)
    opt0     <- array(0, numSim)
    smy0     <- array(0, numSim)
    modified.smy0     <- array(0, numSim)
    
    vsh0     <- array(0, numSim)
    modified.vsh0     <- array(0, numSim)

    feb0     <- array(0, numSim)
    fes0     <- array(0, numSim)
    Feb0     <- array(0, numSim)
    rebayes0 <- array(0, numSim)
    
    f00      <- array(0, numSim)
    f10      <- array(0, numSim)
    f20      <- array(0, numSim)
    
    
    tsse1.sSq <- array(0, numSim) 
    tsse1.ljs <- tsse1.sSq
    tsse1.opt <- tsse1.sSq
    tsse1.smy <- tsse1.sSq
    tsse1.modified.smy <- tsse1.sSq
    
    tsse1.vsh <- tsse1.sSq
    tsse1.modified.vsh <- tsse1.sSq
    tsse1.feb <- tsse1.sSq
    tsse1.fes <- tsse1.sSq
    tsse1.Feb <- tsse1.sSq
    tsse1.rebayes <- tsse1.sSq

    tsse2.sSq <- array(0, numSim)
    tsse2.ljs <- tsse2.sSq
    tsse2.opt <- tsse2.sSq
    tsse2.smy <- tsse2.sSq
    tsse2.modified.smy <- tsse2.sSq

    tsse2.vsh <- tsse2.sSq
    tsse2.modified.vsh <- tsse2.sSq
    tsse2.feb <- tsse2.sSq
    tsse2.fes <- tsse2.sSq
    tsse2.Feb <- tsse2.sSq
    tsse2.rebayes <- tsse2.sSq

    
    tsse0.sSq <- array(0, numSim)
    tsse0.ljs <- tsse0.sSq
    tsse0.opt <- tsse0.sSq
    tsse0.smy <- tsse0.sSq
    tsse0.modified.smy <- tsse0.sSq
    
    tsse0.vsh <- tsse0.sSq
    tsse0.modified.vsh <- tsse0.sSq
    tsse0.feb <- tsse0.sSq
    tsse0.fes <- tsse0.sSq
    tsse0.Feb <- tsse0.sSq
    tsse0.rebayes <- tsse0.sSq
    
    tssep.sSq <- array(0, numSim)
    tssep.ljs <- tssep.sSq
    tssep.opt <- tssep.sSq
    tssep.smy <- tssep.sSq
    tssep.modified.smy <- tssep.sSq
    
    tssep.vsh <- tssep.sSq
    tssep.modified.vsh <- tssep.sSq
    tssep.feb <- tssep.sSq
    tssep.fes <- tssep.sSq
    tssep.Feb <- tssep.sSq
    tssep.rebayes <- tssep.sSq
    
    
    ##-- Simulation ---------------------------------------------------------------------------------------------------------------------------------------------------------
    
    for(i in  1:numSim){
        if( i%% 10==0 ){
          print(date())
          print( paste( "number of simulation:", i, sep="") )
        }
    
        sigmaSq <- GenerateVar(G, prior, a, b)
        sSq        <- sigmaSq*rchisq(G,df)/df

        if(bw=='ucv')
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
        modified.smy        <- modified.SMY(sSq,df)
        
        vsh        <- VSH(sSq,df)
        modified.vsh        <- modified.VSH(sSq,df)

        opt        <- OPT(sSq,df)
        feb        <- fEB(sSq,df,bandwidth)
        fes        <- fES(sSq,df,bandwidth)
        Feb        <- FEB(sSq,df)
        rebayes <- GVmix(sSq, array(df,G) )$dy
        ## feb <- Feb
        ## fes <- Feb
        ## rebayes <- Feb

        
        sigmaSq0 <- GenerateVar(1, prior, a, b)
        sSq0[i]    <- sigmaSq0*rchisq(1,df)/df
        
        sSq_ext    <- c(sSq,sSq0[i])
        
        ## LJS
        B0 <-LJST(sSq,df)[1]; a00 <- LJST(sSq,df)[2]
        
        ## TW    
        OPTT_vec   <- OPTT(sSq,df)
        sSqpool0   <- OPTT_vec[1]
        ahat0      <- OPTT_vec[2]
        nsSqpool0  <- OPTT_vec[3]
        nahat0     <- OPTT_vec[4]
        opt0_temp  <- ((hf(df,1,G)*sSqpool0 )^ahat0 )*(hf(df,1,1)*sSq0[i]^(1-ahat0))
        
    
        ## fEBV, fEBVS
        f00[i]   <- dkde(sSq,sSq0[i],deriv.order=0, h=bandwidth[1])[[8]]
        f10[i]   <- dkde(sSq,sSq0[i],deriv.order=1, h=bandwidth[2])[[8]]
        f20[i]   <- dkde(sSq,sSq0[i],deriv.order=2, h=bandwidth[3])[[8]]
        
        ## Smyth
        s0    <- squeezeVar(sSq, df)$var.prior
        d0    <- squeezeVar(sSq, df)$df.prior
        if (is.infinite(d0)) d0=1000000
        
        ## Estimator for New obs.
        
        ljs0[i]  <- B0*(sSq0[i]^(1-a00))*(sSq0[i]^a00)
        opt0[i]  <- ((hf(df,1,G)*nsSqpool0)^nahat0)*(hf(df,1,1)*opt0_temp^(1-nahat0))
        smy0[i]  <- (d0*s0+df*sSq0[i])/(d0+df)
        
        
        ## Modified Smyth method
        res    <-squeezeVar(sSq, df)
        a0 = res$df.prior/2 # prior parameters
        b0 = (res$df.prior/2)*res$var.prior
        a1 = a0 + df/2
        b1 = b0 + df*sSq0[i]/2
        modified.smy0[i] <- b1/(a1-2) # instead of usual estimate, b1/a1, which is inverse of posterior mean of precision
        
        
        vsh0[i]  <- VSH(sSq_ext,df)[(G+1)]

        vsh <- vash( sqrt(sSq),df )
        a0 = vsh$fitted.g$alpha
        b0 = vsh$fitted.g$beta
        
        pi0 = vsh$fitted.g$pi
        a1 = a0 + df/2
        b1 = b0 + df* sSq0[i]/2
        
        modified.vsh0[i] <-  sum( pi0 * b1/(a1-2 ) )

        
        feb0[i]  <- ifelse(((df-2)/(df*sSq0[i])-2/df*(f10[i]/max(f00[i],0.00001)))^{-1}<0,sSq0[i],((df-2)/(df*sSq0[i])-2/df*(f10[i]/max(f00[i],0.00001)))^{-1})
        Feb0[i]  <- FEB(sSq_ext,df)[(G+1)]
        rebayes0[i] <- rebayes[ which( abs(sSq-sSq0[i])==min( abs(sSq-sSq0[i]) ) ) ]
        
        if (is.na(fes0[i])){
          fes0[i] <- sSq0[i]
        }
        else if (fes0[i]<=0){
          fes0[i] <- sSq0[i]
        }
        else if (abs(fes0[i])<0.0001){
          fes0[i] <- 0.0001
        }
        else if (abs(fes0[i])==Inf){
          fes0[i] <- 10000
        }else {
          fes0[i]  <- (df*(df-2)*sSq0[i]*f00[i]-2*df*sSq0[i]^2*f10[i])/(4*sSq0[i]^2*f20[i]-4*(df-2)*sSq0[i]*f10[i]+df*(df-2)*f00[i])
        }
        
        
        ##----------------------------------------------------------------------------------------------------    
        ## L0
        tsse0.sSq[i] <-(mean((sSq0[i]-sigmaSq0)^2,na.rm=T))  
        tsse0.ljs[i] <-(mean((ljs0[i]-sigmaSq0)^2,na.rm=T))  
        tsse0.opt[i] <-(mean((opt0[i]-sigmaSq0)^2,na.rm=T))  
        tsse0.smy[i] <-(mean((smy0[i]-sigmaSq0)^2,na.rm=T))  
        tsse0.modified.smy[i] <-(mean((modified.smy0[i]-sigmaSq0)^2,na.rm=T))  
        tsse0.vsh[i] <-(mean((vsh0[i]-sigmaSq0)^2,na.rm=T))
        tsse0.modified.vsh[i] <-(mean((modified.vsh0[i]-sigmaSq0)^2,na.rm=T))  
        
        tsse0.feb[i] <-(mean((feb0[i]-sigmaSq0)^2,na.rm=T))  
        tsse0.fes[i] <-(mean((fes0[i]-sigmaSq0)^2,na.rm=T))  
        tsse0.Feb[i] <-(mean((Feb0[i]-sigmaSq0)^2,na.rm=T))
        tsse0.rebayes[i] <-(mean((rebayes0[i]-sigmaSq0)^2,na.rm=T))  
        
    
        ##L1p
        tssep.sSq[i] <-(mean((sSq0[i]/sigmaSq0-1)^2,na.rm=T))
        tssep.ljs[i] <-(mean((ljs0[i]/sigmaSq0-1)^2,na.rm=T))
        tssep.opt[i] <-(mean((opt0[i]/sigmaSq0-1)^2,na.rm=T))
        tssep.smy[i] <-(mean((smy0[i]/sigmaSq0-1)^2,na.rm=T))
        tssep.modified.smy[i] <-(mean((modified.smy0[i]/sigmaSq0-1)^2,na.rm=T))
        
        tssep.vsh[i] <-(mean((vsh0[i]/sigmaSq0-1)^2,na.rm=T))
        tssep.modified.vsh[i] <-(mean((modified.vsh0[i]/sigmaSq0-1)^2,na.rm=T))
        tssep.feb[i] <-(mean((feb0[i]/sigmaSq0-1)^2,na.rm=T))
        tssep.fes[i] <-(mean((fes0[i]/sigmaSq0-1)^2,na.rm=T))
        tssep.Feb[i] <-(mean((Feb0[i]/sigmaSq0-1)^2,na.rm=T))
        tssep.rebayes[i] <-(mean((rebayes0[i]/sigmaSq0-1)^2,na.rm=T))
        
        ##L1
        tsse1.sSq[i] <-(mean((sigmaSq0/sSq0[i]-1)^2,na.rm=T))
        tsse1.ljs[i] <-(mean((sigmaSq0/ljs0[i]-1)^2,na.rm=T))
        tsse1.opt[i] <-(mean((sigmaSq0/opt0[i]-1)^2,na.rm=T))
        tsse1.smy[i] <-(mean((sigmaSq0/smy0[i]-1)^2,na.rm=T))
        tsse1.modified.smy[i] <-(mean((sigmaSq0/modified.smy0[i]-1)^2,na.rm=T))
        
        tsse1.vsh[i] <-(mean((sigmaSq0/vsh0[i]-1)^2,na.rm=T))
        tsse1.modified.vsh[i] <-(mean((sigmaSq0/modified.vsh0[i]-1)^2,na.rm=T))
        tsse1.feb[i] <-(mean((sigmaSq0/feb0[i]-1)^2,na.rm=T))
        tsse1.fes[i] <-(mean((sigmaSq0/fes0[i]-1)^2,na.rm=T))
        tsse1.Feb[i] <-(mean((sigmaSq0/Feb0[i]-1)^2,na.rm=T))
        tsse1.rebayes[i] <-(mean((sigmaSq0/rebayes0[i]-1)^2,na.rm=T))
    
        ##L2
        tsse2.sSq[i] <-(mean((sSq0[i]/sigmaSq0-log(sSq0[i]/sigmaSq0)-1),na.rm=T))
        tsse2.ljs[i] <-(mean((ljs0[i]/sigmaSq0-log(ljs0[i]/sigmaSq0)-1),na.rm=T))
        tsse2.opt[i] <-(mean((opt0[i]/sigmaSq0-log(opt0[i]/sigmaSq0)-1),na.rm=T))
        tsse2.smy[i] <-(mean((smy0[i]/sigmaSq0-log(smy0[i]/sigmaSq0)-1),na.rm=T))
        tsse2.modified.smy[i] <-(mean((modified.smy0[i]/sigmaSq0-log(modified.smy0[i]/sigmaSq0)-1),na.rm=T))
        
        tsse2.vsh[i] <-(mean((vsh0[i]/sigmaSq0-log(vsh0[i]/sigmaSq0)-1),na.rm=T))
        tsse2.modified.vsh[i] <-(mean((modified.vsh0[i]/sigmaSq0-log(modified.vsh0[i]/sigmaSq0)-1),na.rm=T))
        tsse2.feb[i] <-(mean((feb0[i]/sigmaSq0-log(feb0[i]/sigmaSq0)-1),na.rm=T))
        tsse2.fes[i] <-(mean((fes0[i]/sigmaSq0-log(fes0[i]/sigmaSq0)-1),na.rm=T))
        tsse2.Feb[i] <-(mean((Feb0[i]/sigmaSq0-log(Feb0[i]/sigmaSq0)-1),na.rm=T))
        tsse2.rebayes[i] <-(mean((rebayes0[i]/sigmaSq0-log(rebayes0[i]/sigmaSq0)-1),na.rm=T))
      
    }
    
    tmse0.sSq<-mean(tsse0.sSq);  tmse1.sSq<-mean(tsse1.sSq); tmse2.sSq<-mean(tsse2.sSq);  tmsep.sSq<-mean(tssep.sSq);
    tmse0.ljs<-mean(tsse0.ljs);  tmse1.ljs<-mean(tsse1.ljs); tmse2.ljs<-mean(tsse2.ljs);  tmsep.ljs<-mean(tssep.ljs); 
    tmse0.opt<-mean(tsse0.opt);  tmse1.opt<-mean(tsse1.opt); tmse2.opt<-mean(tsse2.opt);  tmsep.opt<-mean(tssep.opt); 
    tmse0.smy<-mean(tsse0.smy);  tmse1.smy<-mean(tsse1.smy); tmse2.smy<-mean(tsse2.smy);  tmsep.smy<-mean(tssep.smy); 
    tmse0.modified.smy<-mean(tsse0.modified.smy);  tmse1.modified.smy<-mean(tsse1.modified.smy); tmse2.modified.smy<-mean(tsse2.modified.smy);  tmsep.modified.smy<-mean(tssep.modified.smy); 
    
    tmse0.vsh<-mean(tsse0.vsh);  tmse1.vsh<-mean(tsse1.vsh); tmse2.vsh<-mean(tsse2.vsh);  tmsep.vsh<-mean(tssep.vsh);
    tmse0.modified.vsh<-mean(tsse0.modified.vsh);  tmse1.modified.vsh<-mean(tsse1.modified.vsh); tmse2.modified.vsh<-mean(tsse2.modified.vsh);  tmsep.modified.vsh<-mean(tssep.modified.vsh); 
    
    tmse0.feb<-mean(tsse0.feb);  tmse1.feb<-mean(tsse1.feb); tmse2.feb<-mean(tsse2.feb);  tmsep.feb<-mean(tssep.feb); 
    tmse0.fes<-mean(tsse0.fes);  tmse1.fes<-mean(tsse1.fes); tmse2.fes<-mean(tsse2.fes);  tmsep.fes<-mean(tssep.fes); 
    tmse0.Feb<-mean(tsse0.Feb);  tmse1.Feb<-mean(tsse1.Feb); tmse2.Feb<-mean(tsse2.Feb);  tmsep.Feb<-mean(tssep.Feb);
    tmse0.rebayes<-mean(tsse0.rebayes);  tmse1.rebayes<-mean(tsse1.rebayes); tmse2.rebayes<-mean(tsse2.rebayes);  tmsep.rebayes<-mean(tssep.rebayes); 
    
    
    tmsd0.sSq<-sd(tsse0.sSq);  tmsd1.sSq<-sd(tsse1.sSq); tmsd2.sSq<-sd(tsse2.sSq);  tmsdp.sSq<-sd(tssep.sSq);
    tmsd0.ljs<-sd(tsse0.ljs);  tmsd1.ljs<-sd(tsse1.ljs); tmsd2.ljs<-sd(tsse2.ljs);  tmsdp.ljs<-sd(tssep.ljs); 
    tmsd0.opt<-sd(tsse0.opt);  tmsd1.opt<-sd(tsse1.opt); tmsd2.opt<-sd(tsse2.opt);  tmsdp.opt<-sd(tssep.opt); 
    tmsd0.smy<-sd(tsse0.smy);  tmsd1.smy<-sd(tsse1.smy); tmsd2.smy<-sd(tsse2.smy);  tmsdp.smy<-sd(tssep.smy); 
    tmsd0.modified.smy<-sd(tsse0.modified.smy);  tmsd1.modified.smy<-sd(tsse1.modified.smy); tmsd2.modified.smy<-sd(tsse2.modified.smy);  tmsdp.modified.smy<-sd(tssep.modified.smy); 
    
    tmsd0.vsh<-sd(tsse0.vsh);  tmsd1.vsh<-sd(tsse1.vsh); tmsd2.vsh<-sd(tsse2.vsh);  tmsdp.vsh<-sd(tssep.vsh);
    tmsd0.modified.vsh<-sd(tsse0.modified.vsh);  tmsd1.modified.vsh<-sd(tsse1.modified.vsh); tmsd2.modified.vsh<-sd(tsse2.modified.vsh);  tmsdp.modified.vsh<-sd(tssep.modified.vsh); 
    tmsd0.feb<-sd(tsse0.feb);  tmsd1.feb<-sd(tsse1.feb); tmsd2.feb<-sd(tsse2.feb);  tmsdp.feb<-sd(tssep.feb); 
    tmsd0.fes<-sd(tsse0.fes);  tmsd1.fes<-sd(tsse1.fes); tmsd2.fes<-sd(tsse2.fes);  tmsdp.fes<-sd(tssep.fes); 
    tmsd0.Feb<-sd(tsse0.Feb);  tmsd1.Feb<-sd(tsse1.Feb); tmsd2.Feb<-sd(tsse2.Feb);  tmsdp.Feb<-sd(tssep.Feb); 
    tmsd0.rebayes<-sd(tsse0.rebayes);  tmsd1.rebayes<-sd(tsse1.rebayes); tmsd2.rebayes<-sd(tsse2.rebayes);  tmsdp.rebayes<-sd(tssep.rebayes); 
    
    
    
    L0  <-round(c(tmse0.sSq, tmse0.ljs, tmse0.opt, tmse0.smy, tmse0.modified.smy, tmse0.vsh, tmse0.modified.vsh, tmse0.feb, tmse0.fes, tmse0.rebayes, tmse0.Feb),3)
    L1  <-round(c(tmse1.sSq, tmse1.ljs, tmse1.opt, tmse1.smy, tmse1.modified.smy, tmse1.vsh, tmse1.modified.vsh, tmse1.feb, tmse1.fes, tmse1.rebayes, tmse1.Feb),3)
    L2  <-round(c(tmse2.sSq, tmse2.ljs, tmse2.opt, tmse2.smy, tmse2.modified.smy, tmse2.vsh, tmse2.modified.vsh, tmse2.feb, tmse2.fes, tmse2.rebayes, tmse2.Feb),3)
    Lp  <-round(c(tmsep.sSq, tmsep.ljs, tmsep.opt, tmsep.smy, tmsep.modified.smy, tmsep.vsh, tmsep.modified.vsh, tmsep.feb, tmsep.fes, tmsep.rebayes, tmsep.Feb),3)
    
    Lsd0<-round(c(tmsd0.sSq, tmsd0.ljs, tmsd0.opt, tmsd0.smy, tmsd0.modified.smy, tmsd0.vsh, tmsd0.modified.vsh, tmsd0.feb,  tmse0.fes, tmsd0.rebayes, tmsd0.Feb),3)
    Lsd1<-round(c(tmsd1.sSq, tmsd1.ljs, tmsd1.opt, tmsd1.smy, tmsd1.modified.smy, tmsd1.vsh, tmsd1.modified.vsh, tmsd1.feb,  tmsd1.fes, tmsd1.rebayes, tmsd1.Feb),3)
    Lsd2<-round(c(tmsd2.sSq, tmsd2.ljs, tmsd2.opt, tmsd2.smy, tmsd2.modified.smy, tmsd2.vsh, tmsd2.modified.vsh, tmsd2.feb,  tmsd2.fes, tmsd2.rebayes, tmsd2.Feb),3)
    Lsdp<-round(c(tmsdp.sSq, tmsdp.ljs, tmsdp.opt, tmsdp.smy, tmsdp.modified.smy, tmsdp.vsh, tmsdp.modified.vsh, tmsdp.feb,  tmsdp.fes, tmsdp.rebayes, tmsdp.Feb),3)
    
    L   <- rbind(L0, L1, L2, Lp)
    Lsd <- rbind(Lsd0, Lsd1, Lsd2, Lsdp)
    
    ##-- Saving RData -------------------------------------------------------------------------------------------------------------------
    res <- list( tsse0.sSq = tsse0.sSq, tsse0.ljs=tsse0.ljs, tsse0.opt=tsse0.opt, tsse0.smy = tsse0.smy, tsse0.modified.smy=tsse0.modified.smy, tsse0.vsh=tsse0.vsh, tsse0.modified.vsh=tsse0.modified.vsh, tsse0.feb=tsse0.feb, tsse0.fes=tsse0.fes, tsse0.Feb=tsse0.Feb, tsse0.rebayes= tsse0.rebayes,  tsse1.sSq = tsse1.sSq, tsse1.ljs=tsse1.ljs, tsse1.opt=tsse1.opt, tsse1.smy = tsse1.smy, tsse1.modified.smy=tsse1.modified.smy, tsse1.vsh=tsse1.vsh, tsse1.modified.vsh=tsse1.modified.vsh, tsse1.feb=tsse1.feb, tsse1.fes=tsse1.fes, tsse1.Feb=tsse1.Feb, tsse1.rebayes= tsse1.rebayes,  tsse2.sSq = tsse2.sSq, tsse2.ljs=tsse2.ljs, tsse2.opt=tsse2.opt, tsse2.smy = tsse2.smy, tsse2.modified.smy=tsse2.modified.smy, tsse2.vsh=tsse2.vsh, tsse2.modified.vsh=tsse2.modified.vsh, tsse2.feb=tsse2.feb, tsse2.fes=tsse2.fes, tsse2.Feb=tsse2.Feb, tsse2.rebayes= tsse2.rebayes,  tssep.sSq = tssep.sSq, tssep.ljs=tssep.ljs, tssep.opt=tssep.opt, tssep.smy = tssep.smy, tssep.modified.smy=tssep.modified.smy, tssep.vsh=tssep.vsh, tssep.modified.vsh=tssep.modified.vsh, tssep.feb=tssep.feb, tssep.fes=tssep.fes, tssep.Feb=tssep.Feb, tssep.rebayes= tssep.rebayes, L0=L0, L1=L1, L2=L2, Lp=Lp, Lsd0=Lsd0, Lsd1=Lsd1, Lsd2=Lsd2, Lsdp=Lsdp )

    path      <- paste("./Data/")
    path1     <- paste("./Figure/")
    
    
    settname  <- paste("FinBay_L1_prior_",prior,"_G_",G,"_a_",a,"_b_",b,sep="")
    filename  <- paste(path, settname, ".Rdata", sep="")
    save(res, file = filename)
    
    
#------------------------------------------------------------------------------------------------------------#
    
  }
  
  if(PLOT==TRUE)
    {
      path      <- paste("./Data/")
      path1     <- paste("./Figure/")
      
      
      settname  <- paste("FinBay_L1_prior_",prior,"_G_",G,"_a_",a,"_b_",b,sep="")
      filename  <- paste(path, settname, ".Rdata", sep="")
      load(filename)
      
      
    }

  path1     <- paste("./Figure/")
  ##-- Saving Figures ----------------------------------------------------------------------------------------------------------------------------
  
  if (prior==1){
    PRIOR<-"Inverse Gamma"
  } else if (prior==2){
    PRIOR<-"Lognormal"
  } else if (prior==3){
    PRIOR<-"Mix InGamma"
  } else if (prior==4){
    PRIOR<-"Mix Lognormal"
  } else if (prior==5){
    PRIOR <- "InvGauss"
  } else if (prior==6){
    PRIOR <- "Mix InvGauss"
  } else if (prior==7){
    PRIOR <- "Two Point"
  }

  
  namelist=c("sSq",	"ELJS",	"TW",	"Smyth", "Modified Smyth",	"Vash",	"Modified Vash", "fEBV",	"fEBVS", "REBayes", "FEBV")
  
  
  figurename1 <- paste(path1,"FinBay_",settname,"_L1.pdf",sep="")
  pdf(figurename1,width=8,height=5)
  barplot(log(L1,10), ylab="log(MSE)", main=bquote(paste("Finite Byaes:" ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
  dev.off()
  
  
  
                                        ##-- Showing Figures ----------------------------------------------------------------------------------------------------------------------------
  
  barplot(log(L1,10), ylab="MSE (log_scale)", main=bquote(paste("[L1, 1%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
  

  
}
