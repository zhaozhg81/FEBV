
#-- Saving Figures ----------------------------------------------------------------------------------------------------------------------------

if (prior==1){
  PRIOR<-"Inverse Gamma"
} else if (prior==2){
  PRIOR<-"Lognormal"
} else if (prior==3){
  PRIOR<-"Gamma"
} else if (prior==4){
  PRIOR<-"IG Mixture"
} else if (prior==5){
  PRIOR<-"Gamma Mixture"
}

par(mar=c(4, 4, 3, 8))
xloc=-0.2

if (prior<4){
# Comparison of Estimators for theta
  figurename1 <- paste(path1,"MSE_mean_comp_",settname,".pdf",sep="")
  pdf(figurename1,width=8,height=5)
  par(mar=c(5,5,1,10) )
  plot( Result.unknown[,1],Result.unknown[,2],xlab="M",ylab="MSE", ylim=c(0,max(Result.unknown[,2])), xlim=c(0,1),col="white",main=bquote(paste(.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df)," (top ",.(perc*100),"%)")), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( Result.unknown[,1],Result.unknown[,2], pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( Result.unknown[,1],Result.unknown[,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result.unknown[,1],Result.unknown[,4], pch=7 , type="b", col="purple") # TW
  lines( Result.unknown[,1],Result.unknown[,5], pch=17, type="b", col="green4") # Smyth
  lines( Result.unknown[,1],Result.unknown[,6], pch=16, type="b", col="gold4" ) # Vash
  lines( Result.unknown[,1],Result.unknown[,7], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result.unknown[,1],Result.unknown[,8], pch=8 , type="b", col="orange") # fEBVS
  lines( Result.unknown[,1],Result.unknown[,9], pch=19, type="b", col="red"   ) # FEBV
  par(xpd=TRUE)
  legend("topright",horiz=F, inset = c(-0.3,0.2), legend=c("sSq","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","F-EBV"),col=c("cyan3","blue","purple","green4","gold4", "brown","orange","red"),pch=c(3,15,7,17,16,9,8,19), ncol=1, lty=1:1, cex=1)
  dev.off()
  
  
  #-- Showing Figures ----------------------------------------------------------------------------------------------------------------------------
  
  plot( Result.unknown[,1],Result.unknown[,2],xlab="M",ylab="Coverage Prob", ylim=c(0,max(Result.unknown[,2])), xlim=c(0,1),col="white",main=bquote(paste("[MSE for Selected Mean]   ",.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df)," (top ",.(perc*100),"%)")), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( Result.unknown[,1],Result.unknown[,2], pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( Result.unknown[,1],Result.unknown[,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result.unknown[,1],Result.unknown[,4], pch=7 , type="b", col="purple") # TW
  lines( Result.unknown[,1],Result.unknown[,5], pch=17, type="b", col="green4") # Smyth
  lines( Result.unknown[,1],Result.unknown[,6], pch=16, type="b", col="gold4" ) # Vash
  lines( Result.unknown[,1],Result.unknown[,7], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result.unknown[,1],Result.unknown[,8], pch=8 , type="b", col="orange") # fEBVS
  lines( Result.unknown[,1],Result.unknown[,9], pch=19, type="b", col="red"   ) # FEBV
  legend("topright",horiz=F, inset = c(xloc,0), legend=c("sSq","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","F-EBV"),col=c("cyan3","blue","purple","green4","gold4", "brown","orange","red"),pch=c(3,15,7,17,16,9,8,19), ncol=1, lty=1:1, cex=1)
}





if (prior>3){
  # Comparison of Estimators for theta
  figurename1 <- paste(path1,"MSE_mean_comp_",settname,".pdf",sep="")
  pdf(figurename1,width=8,height=5)
  par(mar=c(5,5,1,10) )
  plot( Result.unknown[,1],Result.unknown[,2],xlab="M",ylab="MSE", ylim=c(0,max(Result.unknown[,2])), xlim=c(0,1),col="white",main=bquote(paste("[MSE for Selected Mean]   ",.(PRIOR),"    a=",.(a),",  b=",.(b),",  df=",.(df)," (top ",.(perc*100),"%)")), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( Result.unknown[,1],Result.unknown[,2], pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( Result.unknown[,1],Result.unknown[,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result.unknown[,1],Result.unknown[,4], pch=7 , type="b", col="purple") # TW
  lines( Result.unknown[,1],Result.unknown[,5], pch=17, type="b", col="green4") # Smyth
  lines( Result.unknown[,1],Result.unknown[,6], pch=16, type="b", col="gold4" ) # Vash
  lines( Result.unknown[,1],Result.unknown[,7], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result.unknown[,1],Result.unknown[,8], pch=8 , type="b", col="orange") # fEBVS
  lines( Result.unknown[,1],Result.unknown[,9], pch=19, type="b", col="red"   ) # FEBV
  ## legend("topright",horiz=F, inset = c(xloc,0), legend=c("sSq","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","F-EBV"),col=c("cyan3","blue","purple","green4","gold4", "brown","orange","red"),pch=c(3,15,7,17,16,9,8,19), ncol=1, lty=1:1, cex=1)
  dev.off()
  
  
  #-- Showing Figures ----------------------------------------------------------------------------------------------------------------------------
  
  plot( Result.unknown[,1],Result.unknown[,2],xlab="M",ylab="Coverage Prob", ylim=c(0,max(Result.unknown[,2])), xlim=c(0,1),col="white",main=bquote(paste("[MSE for Selected Mean]   ",.(PRIOR),"    a=",.(a),",  b=",.(b),",  df=",.(df)," (top ",.(perc*100),"%)")), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( Result.unknown[,1],Result.unknown[,2], pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( Result.unknown[,1],Result.unknown[,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result.unknown[,1],Result.unknown[,4], pch=7 , type="b", col="purple") # TW
  lines( Result.unknown[,1],Result.unknown[,5], pch=17, type="b", col="green4") # Smyth
  lines( Result.unknown[,1],Result.unknown[,6], pch=16, type="b", col="gold4" ) # Vash
  lines( Result.unknown[,1],Result.unknown[,7], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result.unknown[,1],Result.unknown[,8], pch=8 , type="b", col="orange") # fEBVS
  lines( Result.unknown[,1],Result.unknown[,9], pch=19, type="b", col="red"   ) # FEBV
  legend("topright",horiz=F, inset = c(xloc,0), legend=c("sSq","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","F-EBV"),col=c("cyan3","blue","purple","green4","gold4", "brown","orange","red"),pch=c(3,15,7,17,16,9,8,19), ncol=1, lty=1:1, cex=1)
}























