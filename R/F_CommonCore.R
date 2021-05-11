
#-- Common Core ----------------------------------------------------------------------------------------------------------------------------------------------
     df      <- n-1   # Degree of Freedom
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
     
     CI.bon.sel.unknown       <- array(0, c(numSim,Ncol))
     CI.sSq.sel.unknown       <- array(0, c(numSim,Ncol))
     CI.ljs.sel.unknown       <- array(0, c(numSim,Ncol))
     CI.opt.sel.unknown       <- array(0, c(numSim,Ncol))
     CI.smy.sel.unknown       <- array(0, c(numSim,Ncol))
     CI.vsh.sel.unknown       <- array(0, c(numSim,Ncol))
     CI.feb.sel.unknown       <- array(0, c(numSim,Ncol))
     CI.fes.sel.unknown       <- array(0, c(numSim,Ncol))
     CI.reb.sel.unknown       <- array(0, c(numSim,Ncol))
     CI.Feb.sel.unknown       <- array(0, c(numSim,Ncol))
     
      tauSqhat.sSq                 <- array(0, numSim)
      tauSqhat.ljs                 <- array(0, numSim)
      tauSqhat.opt                 <- array(0, numSim)
      tauSqhat.smy                 <- array(0, numSim)
      tauSqhat.vsh                 <- array(0, numSim)
      tauSqhat.feb                 <- array(0, numSim)
      tauSqhat.fes                 <- array(0, numSim)
      tauSqhat.reb                 <- array(0, numSim)
      tauSqhat.Feb                 <- array(0, numSim)

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

     Result.unknown           <- array(0, c(N,41))
     Result.avg.unknown       <- array(0, c(N,41))
     
     M=seq(0,1,0.999/N)
     M=M[2:(N+1)]
     
     
#-- Simulation ----------------------------------------------------------------------------------------------------------------------------------------------

    for(t in  1:N){
        print( paste("t=",t,sep=""))
        for(i in  1:numSim){
           if( i%% 50==0 ){
             print(date())
             print( paste( "number of simulation:", i, sep="") )
           }

          sigmaSq <- GenerateVar(G, 1, 10, 1)
          

          tauSq      <- M[t]/(1-M[t])
          theta      <- rnorm(G,mu,sqrt(tauSq))

           x          <- matrix(rnorm(n*G,0,1),G,n)
##           L          <- matrix(rep(sqrt(sigmaSq),G),G,G)*B
           L          <- diag( sqrt(sigmaSq) ) %*% B
           
           temp       <- L %*% x + array(theta,c(G,n))
           xi         <- apply(temp,1,mean)
           sSq        <- apply(temp,1,var)

           muhat      <- mean(xi) #sum((xi/sSq))/sum(1/sSq)
           T.stat     <- xi/sqrt(sSq)
           max.ind    <- order(abs(T.stat),decreasing=TRUE)[1:1]

           if( i==1 )
             {
               bandwidth <- array(0, 3)
               bandwidth[1] <- h.ucv(sSq, deriv.order=0)$h
               bandwidth[2] <- h.ucv(sSq, deriv.order=1)$h
               bandwidth[3] <- h.ucv(sSq, deriv.order=2)$h                
             }

           bon        <- sSq/n
           sSq        <- sSq/n
           ljs        <- LJS(sSq,df)/n
           opt        <- OPT(sSq,df)/n
           smy        <- SMY(sSq,df)/n
           #vsh        <- VSH(sSq,df)/n
           #feb        <- fEB(sSq,df,bandwidth)/n
           #fes        <- fES(sSq,df,bandwidth)/n
          ## reb        <- GVmix(sSq,array(df,G))$dy/n
           Feb        <- FEB(sSq,df)/n
           reb = Feb
           vsh=Feb
           feb=Feb
           fes=Feb
           
           tauSqhat[i]<- tauSqfn(sSq/n,xi,muhat,alpha)

           tauSqhat.sSq[i]<- tauSqhat[i]
          tauSqhat.ljs[i]<- tauSqhat[i]
          tauSqhat.opt[i]<- tauSqhat[i]
          tauSqhat.smy[i]<- tauSqhat[i]
          tauSqhat.vsh[i]<- tauSqhat[i]
          tauSqhat.feb[i]<- tauSqhat[i]
          tauSqhat.fes[i]<- tauSqhat[i]
          tauSqhat.reb[i]<- tauSqhat[i]
          tauSqhat.Feb[i]<- tauSqhat[i]
          

           CI.bon.unknown   <- CI.bonferroni(sSq/n,df,xi,theta,sigmaSq,muhat,alpha)
           CI.sSq.unknown   <- CI.t(sSq/n,df,xi,theta,sigmaSq,muhat,alpha)
           CI.ljs.unknown   <- CI.est(ljs,tauSqhat.ljs[i],xi,theta,sigmaSq,muhat,alpha)
           CI.opt.unknown   <- CI.est(opt,tauSqhat.opt[i],xi,theta,sigmaSq,muhat,alpha)
           CI.smy.unknown   <- CI.est(smy,tauSqhat.smy[i],xi,theta,sigmaSq,muhat,alpha)
           CI.vsh.unknown   <- CI.est(vsh,tauSqhat.vsh[i],xi,theta,sigmaSq,muhat,alpha)
           CI.feb.unknown   <- CI.est(feb,tauSqhat.feb[i],xi,theta,sigmaSq,muhat,alpha)
           CI.fes.unknown   <- CI.est(fes,tauSqhat.fes[i],xi,theta,sigmaSq,muhat,alpha)
           CI.reb.unknown   <- CI.est(reb,tauSqhat.reb[i],xi,theta,sigmaSq,muhat,alpha)
           CI.Feb.unknown   <- CI.est(Feb,tauSqhat.Feb[i],xi,theta,sigmaSq,muhat,alpha)
           
           # CI : tauSq unknown - selected
           CI.bon.sel.unknown[i,]   <- CI.bon.unknown[max.ind,]
           CI.sSq.sel.unknown[i,]   <- CI.sSq.unknown[max.ind,]
           CI.ljs.sel.unknown[i,]   <- CI.ljs.unknown[max.ind,]
           CI.opt.sel.unknown[i,]   <- CI.opt.unknown[max.ind,]
           CI.smy.sel.unknown[i,]   <- CI.smy.unknown[max.ind,]
           CI.vsh.sel.unknown[i,]   <- CI.vsh.unknown[max.ind,]
           CI.feb.sel.unknown[i,]   <- CI.feb.unknown[max.ind,]
           CI.fes.sel.unknown[i,]   <- CI.fes.unknown[max.ind,]
           CI.reb.sel.unknown[i,]   <- CI.reb.unknown[max.ind,]
           CI.Feb.sel.unknown[i,]   <- CI.Feb.unknown[max.ind,]
        }

      CI.bon.sel.unknown   <-changenames(CI.bon.sel.unknown)
      CI.sSq.sel.unknown   <-changenames(CI.sSq.sel.unknown)
      CI.ljs.sel.unknown   <-changenames(CI.ljs.sel.unknown)
      CI.opt.sel.unknown   <-changenames(CI.opt.sel.unknown)
      CI.smy.sel.unknown   <-changenames(CI.smy.sel.unknown)
      CI.vsh.sel.unknown   <-changenames(CI.vsh.sel.unknown)
      CI.feb.sel.unknown   <-changenames(CI.feb.sel.unknown)
      CI.fes.sel.unknown   <-changenames(CI.fes.sel.unknown)
      CI.reb.sel.unknown   <-changenames(CI.reb.sel.unknown)
      CI.Feb.sel.unknown   <-changenames(CI.Feb.sel.unknown)
      
      bon.sel.unknown.Result=apply(CI.bon.sel.unknown[,7:10],2,mean);
      sSq.sel.unknown.Result=apply(CI.sSq.sel.unknown[,7:10],2,mean);
      ljs.sel.unknown.Result=apply(CI.ljs.sel.unknown[,7:10],2,mean);
      opt.sel.unknown.Result=apply(CI.opt.sel.unknown[,7:10],2,mean);
      smy.sel.unknown.Result=apply(CI.smy.sel.unknown[,7:10],2,mean);
      vsh.sel.unknown.Result=apply(CI.vsh.sel.unknown[,7:10],2,mean);
      feb.sel.unknown.Result=apply(CI.feb.sel.unknown[,7:10],2,mean);
      fes.sel.unknown.Result=apply(CI.fes.sel.unknown[,7:10],2,mean);
      reb.sel.unknown.Result=apply(CI.reb.sel.unknown[,7:10],2,mean);
      Feb.sel.unknown.Result=apply(CI.Feb.sel.unknown[,7:10],2,mean);
      
      Result.unknown[t,]=c(M[t],bon=bon.sel.unknown.Result,sSq=sSq.sel.unknown.Result,ljs=ljs.sel.unknown.Result,opt=opt.sel.unknown.Result,smy=smy.sel.unknown.Result,vsh=vsh.sel.unknown.Result,feb=feb.sel.unknown.Result,fes=fes.sel.unknown.Result,reb=reb.sel.unknown.Result,Feb=Feb.sel.unknown.Result)
    }

#-----------------------------------------------------------------------------------------------------------------------------------------------------
