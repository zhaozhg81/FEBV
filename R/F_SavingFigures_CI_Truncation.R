
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

# Coverage 
  figurename1 <- paste(path1,"Trunc_Coverage_",settname,".pdf",sep="")
  pdf(figurename1,width=8,height=5)
  Llim=min(Result.unknown[,2],Result.unknown[,10],Result.unknown[,14],Result.unknown[,18],Result.unknown[,22])-0.1
  plot( Result.unknown[,1],Result.unknown[,2],xlab="M",ylab="Coverage Prob", ylim=c(0,1), xlim=c(0,1),col="white",main=bquote(paste("[Coverage Prob.]   ",.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df) )), cex.lab=1.3, cex.axis=1.3, cex.main=1.5)
  lines( Result.unknown[,1],rep(0.95-qnorm(0.975)*sqrt(0.05*0.95/1000),N), type="b", col="black",lwd=3,lty=2)
  lines( Result.unknown[,1],Result.unknown[,2] , pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( Result.unknown[,1],Result.unknown[,10], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result.unknown[,1],Result.unknown[,14], pch=7 , type="b", col="purple") # TW
  lines( Result.unknown[,1],Result.unknown[,18], pch=17, type="b", col="green4") # Smyth
  lines( Result.unknown[,1],Result.unknown[,22], pch=16, type="b", col="gold4" ) # Vash
  lines( Result.unknown[,1],Result.unknown[,26], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result.unknown[,1],Result.unknown[,30], pch=8 , type="b", col="orange") # fEBVS
  lines( Result.unknown[,1],Result.unknown[,34], pch=19, type="b", col="red"   ) # FEBV
  legend("bottomright",horiz=F, legend=c("Error Bar","Bonferroni","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","F-EBV","SOC"),col=c("black","cyan3","blue","purple","green4","gold4", "brown","orange","red","green"),pch=c(NA,3,15,7,17,16,9,8,19,10), ncol=1, lty=1:1, cex=1.3)
  dev.off()
  
  # Length Ratio
  figurename2 <- paste(path1,"Trunc_LengRatio_",settname,".pdf",sep="")
  pdf(figurename2,width=8,height=5)
  Llim=max(Result.unknown[,11]/Result.unknown[,3],Result.unknown[,15]/Result.unknown[,3],Result.unknown[,35]/Result.unknown[,3])+0.1
  plot( Result.unknown[,1],Result.unknown[,3],xlab="M",ylab="Length Ratio", ylim=c(0,Llim+0.2), xlim=c(0,1),col="white",main=bquote(paste("[Length Ratio]   ",.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df) )), cex.lab=1.3, cex.axis=1.3, cex.main=1.5)
  lines( Result.unknown[,1],Result.unknown[,11]/Result.unknown[,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result.unknown[,1],Result.unknown[,15]/Result.unknown[,3], pch=7 , type="b", col="purple") # TW
  lines( Result.unknown[,1],Result.unknown[,19]/Result.unknown[,3], pch=17, type="b", col="green4") # Smyth
  lines( Result.unknown[,1],Result.unknown[,23]/Result.unknown[,3], pch=16, type="b", col="gold4" ) # Vash
  lines( Result.unknown[,1],Result.unknown[,27]/Result.unknown[,3], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result.unknown[,1],Result.unknown[,31]/Result.unknown[,3], pch=8 , type="b", col="orange") # fEBVS
  lines( Result.unknown[,1],Result.unknown[,35]/Result.unknown[,3], pch=19, type="b", col="red"   ) # FEBV
  legend("bottomright",horiz=F, legend=c("Error Bar","Bonferroni","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","F-EBV","SOC"),col=c("black","cyan3","blue","purple","green4","gold4", "brown","orange","red","green"),pch=c(NA,3,15,7,17,16,9,8,19,10), ncol=1, lty=1:1, cex=1.3)
  dev.off()
  
  #-- Showing Figures ----------------------------------------------------------------------------------------------------------------------------
  
  # Coverage 
  Llim=min(Result.unknown[,2],Result.unknown[,10],Result.unknown[,14],Result.unknown[,18],Result.unknown[,22])-0.1
  plot( Result.unknown[,1],Result.unknown[,2],xlab="M",ylab="Coverage Prob", ylim=c(0,1), xlim=c(0,1),col="white",main=bquote(paste("[Coverage Prob.]   ",.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df) )), cex.lab=1.3, cex.axis=1.3, cex.main=1.5)
  lines( Result.unknown[,1],rep(0.95-qnorm(0.975)*sqrt(0.05*0.95/1000),N), type="l", col="black",lwd=3,lty=2)
  lines( Result.unknown[,1],Result.unknown[,2] , pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( Result.unknown[,1],Result.unknown[,10], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result.unknown[,1],Result.unknown[,14], pch=7 , type="b", col="purple") # TW
  lines( Result.unknown[,1],Result.unknown[,18], pch=17, type="b", col="green4") # Smyth
  lines( Result.unknown[,1],Result.unknown[,22], pch=16, type="b", col="gold4" ) # Vash
  lines( Result.unknown[,1],Result.unknown[,26], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result.unknown[,1],Result.unknown[,30], pch=8 , type="b", col="orange") # fEBVS
  lines( Result.unknown[,1],Result.unknown[,34], pch=19, type="b", col="red"   ) # FEBV
  legend("bottomright",horiz=F, legend=c("Error Bar","Bonferroni","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","F-EBV","SOC"),col=c("black","cyan3","blue","purple","green4","gold4", "brown","orange","red","green"),pch=c(NA,3,15,7,17,16,9,8,19,10), ncol=1, lty=1:1, cex=1.3)
  
  # Length Ratio
  Llim=max(Result.unknown[,11]/Result.unknown[,3],Result.unknown[,15]/Result.unknown[,3],Result.unknown[,35]/Result.unknown[,3])+0.1
  plot( Result.unknown[,1],Result.unknown[,3],xlab="M",ylab="Length Ratio", ylim=c(0,Llim+0.2), xlim=c(0,1),col="white",main=bquote(paste("[Length Ratio]   ",.(PRIOR), "    a=",.(a),",  b=",.(b),",  df=",.(df) )), cex.lab=1.3, cex.axis=1.3, cex.main=1.5)
  lines( Result.unknown[,1],Result.unknown[,11]/Result.unknown[,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result.unknown[,1],Result.unknown[,15]/Result.unknown[,3], pch=7 , type="b", col="purple") # TW
  lines( Result.unknown[,1],Result.unknown[,19]/Result.unknown[,3], pch=17, type="b", col="green4") # Smyth
  lines( Result.unknown[,1],Result.unknown[,23]/Result.unknown[,3], pch=16, type="b", col="gold4" ) # Vash
  lines( Result.unknown[,1],Result.unknown[,27]/Result.unknown[,3], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result.unknown[,1],Result.unknown[,31]/Result.unknown[,3], pch=8 , type="b", col="orange") # fEBVS
  lines( Result.unknown[,1],Result.unknown[,35]/Result.unknown[,3], pch=19, type="b", col="red"   ) # FEBV
  legend("bottomright",horiz=F, legend=c("Error Bar","Bonferroni","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","F-EBV","SOC"),col=c("black","cyan3","blue","purple","green4","gold4", "brown","orange","red","green"),pch=c(NA,3,15,7,17,16,9,8,19,10), ncol=1, lty=1:1, cex=1.3)
  
  
  