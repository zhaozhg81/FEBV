
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

SOC<-read.csv("./ZEMCISel/SOC.txt", header=FALSE)

if (prior==1 && a==2.44  ) {
  socresult<-SOC[,1:2]
} else if (prior==1 && a==3     ){
  socresult<-SOC[,3:4]
} else if (prior==1 && a==6     ){
  socresult<-SOC[,5:6]
} else if (prior==1 && a==18    ){
  socresult<-SOC[,7:8]
} else if (prior==2 && a==-0.723){
  socresult<-SOC[,9:10]
} else if (prior==2 && a==-0.805){
  socresult<-SOC[,11:12]
} else if (prior==2 && a==-1.040){
  socresult<-SOC[,13:14]
} else if (prior==2 && a==-1.282){
  socresult<-SOC[,15:16]
} else if (prior==3 && a==1     ){
  socresult<-SOC[,17:18]
} else if (prior==3 && a==1.78  ){
  socresult<-SOC[,19:20]
} else if (prior==3 && a==4     ){
  socresult<-SOC[,21:22]
} else if (prior==3 && a==16    ){
  socresult<-SOC[,23:24]
} else if (prior==4){
  socresult<-SOC[,25:26]
} else if (prior==5){
  socresult<-SOC[,27:28]
}
SOCR<-socresult
Result.unknown=cbind(Result.unknown,socresult)
par(mar=c(4, 4, 3, 8))

# Coverage 
  figurename1 <- paste(path1,"Coverage_",settname,".pdf",sep="")
  pdf(figurename1,width=8,height=5)
  Llim=min(Result.unknown[,2],Result.unknown[,10],Result.unknown[,14],Result.unknown[,18],Result.unknown[,22])-0.1
  plot( Result.unknown[,1],Result.unknown[,2],xlab="M",ylab="Coverage Prob", ylim=c(0,1), xlim=c(0,1),col="white",main=bquote(paste("[Coverage Prob.]   ",.(PRIOR), ",  df=",.(df) )), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( Result.unknown[,1],rep(0.95-qnorm(0.975)*sqrt(0.05*0.95/1000),N), type="l", col="black",lwd=3,lty=2)
  lines( Result.unknown[,1],Result.unknown[,2] , pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( Result.unknown[,1],Result.unknown[,10], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result.unknown[,1],Result.unknown[,14], pch=7 , type="b", col="purple") # TW
  lines( Result.unknown[,1],Result.unknown[,18], pch=17, type="b", col="green4") # Smyth
  lines( Result.unknown[,1],Result.unknown[,22], pch=16, type="b", col="gold4" ) # Vash
  lines( Result.unknown[,1],Result.unknown[,26], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result.unknown[,1],Result.unknown[,30], pch=8 , type="b", col="orange") # fEBVS
  lines( Result.unknown[,1],Result.unknown[,34], pch=19, type="b", col="red"   ) # FEBV
  lines( Result.unknown[,1],Result.unknown[,38], pch=10, type="b", col="green ") # SOC
  legend("topright",horiz=F, inset = c(xloc,0), legend=c("Error Bar","Bonferroni","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","F-EBV","SOC"),col=c("black","cyan3","blue","purple","green4","gold4", "brown","orange","red","green"),pch=c(NA,3,15,7,17,16,9,8,19,10), ncol=1, lty=1:1, cex=1)
  dev.off()
  
  # Length Ratio
  figurename2 <- paste(path1,"LengRatio_",settname,".pdf",sep="")
  pdf(figurename2,width=8,height=5)
  Llim=max(Result.unknown[,11]/Result.unknown[,3],Result.unknown[,15]/Result.unknown[,3],Result.unknown[,35]/Result.unknown[,3])+0.1
  plot( Result.unknown[,1],Result.unknown[,3],xlab="M",ylab="Length Ratio", ylim=c(0,Llim+0.2), xlim=c(0,1),col="white",main=bquote(paste("[Length Ratio]   ",.(PRIOR),",  df=",.(df) )), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( Result.unknown[,1],Result.unknown[,11]/Result.unknown[,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result.unknown[,1],Result.unknown[,15]/Result.unknown[,3], pch=7 , type="b", col="purple") # TW
  lines( Result.unknown[,1],Result.unknown[,19]/Result.unknown[,3], pch=17, type="b", col="green4") # Smyth
  lines( Result.unknown[,1],Result.unknown[,23]/Result.unknown[,3], pch=16, type="b", col="gold4" ) # Vash
  lines( Result.unknown[,1],Result.unknown[,27]/Result.unknown[,3], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result.unknown[,1],Result.unknown[,31]/Result.unknown[,3], pch=8 , type="b", col="orange") # fEBVS
  lines( Result.unknown[,1],Result.unknown[,35]/Result.unknown[,3], pch=19, type="b", col="red"   ) # FEBV
  lines( Result.unknown[,1],Result.unknown[,39]/Result.unknown[,3], pch=10, type="b", col="green ") # SOC
  legend("topright",horiz=F, inset = c(xloc,0), legend=c("ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","F-EBV","SOC"),col=c("blue","purple","green4","gold4", "brown","orange","red","green"),pch=c(15,7,17,16,9,8,19,10), ncol=1, lty=1:1, cex=1)
  dev.off()
  
  #-- Showing Figures ----------------------------------------------------------------------------------------------------------------------------
  
  # Coverage 
  Llim=min(Result.unknown[,2],Result.unknown[,10],Result.unknown[,14],Result.unknown[,18],Result.unknown[,22])-0.1
  plot( Result.unknown[,1],Result.unknown[,2],xlab="M",ylab="Coverage Prob", ylim=c(0,1), xlim=c(0,1),col="white",main=bquote(paste("[Coverage Prob.]   ",.(PRIOR), ",  df=",.(df) )), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( Result.unknown[,1],rep(0.95-qnorm(0.975)*sqrt(0.05*0.95/1000),N), type="l", col="black",lwd=3,lty=2)
  lines( Result.unknown[,1],Result.unknown[,2] , pch=3 , type="b", col="cyan3" ) # Bonferroni
  lines( Result.unknown[,1],Result.unknown[,10], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result.unknown[,1],Result.unknown[,14], pch=7 , type="b", col="purple") # TW
  lines( Result.unknown[,1],Result.unknown[,18], pch=17, type="b", col="green4") # Smyth
  lines( Result.unknown[,1],Result.unknown[,22], pch=16, type="b", col="gold4" ) # Vash
  lines( Result.unknown[,1],Result.unknown[,26], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result.unknown[,1],Result.unknown[,30], pch=8 , type="b", col="orange") # fEBVS
  lines( Result.unknown[,1],Result.unknown[,34], pch=19, type="b", col="red"   ) # FEBV
  lines( Result.unknown[,1],Result.unknown[,38], pch=10, type="b", col="green ") # SOC
  legend("topright",horiz=F, inset = c(xloc,0), legend=c("Error Bar","Bonferroni","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","F-EBV","SOC"),col=c("black","cyan3","blue","purple","green4","gold4", "brown","orange","red","green"),pch=c(NA,3,15,7,17,16,9,8,19,10), ncol=1, lty=1:1, cex=1)
  
  # Length Ratio
  Llim=max(Result.unknown[,11]/Result.unknown[,3],Result.unknown[,15]/Result.unknown[,3],Result.unknown[,35]/Result.unknown[,3])+0.1
  plot( Result.unknown[,1],Result.unknown[,3],xlab="M",ylab="Length Ratio", ylim=c(0,Llim+0.2), xlim=c(0,1),col="white",main=bquote(paste("[Length Ratio]   ",.(PRIOR),",  df=",.(df) )), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
  lines( Result.unknown[,1],Result.unknown[,11]/Result.unknown[,3], pch=15, type="b", col="blue"  ) # ELJS
  lines( Result.unknown[,1],Result.unknown[,15]/Result.unknown[,3], pch=7 , type="b", col="purple") # TW
  lines( Result.unknown[,1],Result.unknown[,19]/Result.unknown[,3], pch=17, type="b", col="green4") # Smyth
  lines( Result.unknown[,1],Result.unknown[,23]/Result.unknown[,3], pch=16, type="b", col="gold4" ) # Vash
  lines( Result.unknown[,1],Result.unknown[,27]/Result.unknown[,3], pch=9 , type="b", col="brown" ) # fEBV
  lines( Result.unknown[,1],Result.unknown[,31]/Result.unknown[,3], pch=8 , type="b", col="orange") # fEBVS
  lines( Result.unknown[,1],Result.unknown[,35]/Result.unknown[,3], pch=19, type="b", col="red"   ) # FEBV
  lines( Result.unknown[,1],Result.unknown[,39]/Result.unknown[,3], pch=10, type="b", col="green ") # SOC
  legend("topright",horiz=F, inset = c(xloc,0), legend=c("ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","F-EBV","SOC"),col=c("blue","purple","green4","gold4", "brown","orange","red","green"),pch=c(15,7,17,16,9,8,19,10), ncol=1, lty=1:1, cex=1)
  
  