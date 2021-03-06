library(latex2exp)

#-- Figures ----------------------------------------------------------------------------------------------------------------------------------------------

if (Matrix==1){
  MAT<-"AR"
} else if (Matrix==2){
  MAT<-"Band"
} else if (Matrix==3){
  MAT<-"Sparse"
} 

par(mar=c(4, 4, 3, 8))
xloc=-0.3

## Coverage
figurename1 <- paste(path1,"Coverage_",settname,".pdf",sep="")
pdf(figurename1,width=8,height=5)
par(mar=c(5,5,1,10) )
Llim=min(Result.unknown[,2],Result.unknown[,10],Result.unknown[,14],Result.unknown[,18],Result.unknown[,22])-0.1
if(Matrix==1){
  plot( Result.unknown[,1],Result.unknown[,2], xlab=TeX(r'($\tau^2/(1+\tau^2)$ )'),ylab=TeX(r'(Coverage Prob. of $\theta_{(1)}$)'), ylim=c(0,1), xlim=c(0,1),col="white",main=bquote(paste("AR    ", ~rho,"=",.(rho),",   df=",.(df))), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
}
if(Matrix==2)
  plot( Result.unknown[,1],Result.unknown[,2],xlab=TeX(r'($\tau^2/(1+\tau^2)$ )'),ylab=TeX(r'(Coverage Prob. of $\theta_{(1)}$)'), ylim=c(0,1), xlim=c(0,1),col="white",main=bquote(paste("Banded    ", ~rho,"=",.(rho),",   k=",.(k),",   df=",.(df))), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
if(Matrix==3)
  plot( Result.unknown[,1],Result.unknown[,2], xlab=TeX(r'($\tau^2/(1+\tau^2)$ )'),ylab=TeX(r'(Coverage Prob. of $\theta_{(1)}$)'), ylim=c(0,1), xlim=c(0,1),col="white",main=bquote(paste("Sparse    ", ~rho," in (",.(a),",",.(b),"),   df=",.(df),",   Sparsity=",.(Sparsity),"%")), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)

lines( Result.unknown[,1],rep(0.95-qnorm(0.975)*sqrt(0.05*0.95/1000),N), type="b", col="black",lwd=3,lty=2)
lines( Result.unknown[,1],Result.unknown[,2] , pch=3 , type="b", col="cyan3" ) # Bonferroni
lines( Result.unknown[,1],Result.unknown[,10], pch=15, type="b", col="blue"  ) # ELJS
lines( Result.unknown[,1],Result.unknown[,14], pch=7 , type="b", col="purple") # TW
lines( Result.unknown[,1],Result.unknown[,18], pch=17, type="b", col="green4") # Smyth
lines( Result.unknown[,1],Result.unknown[,22], pch=16, type="b", col="gold4" ) # Vash
lines( Result.unknown[,1],Result.unknown[,26], pch=9 , type="b", col="brown" ) # fEBV
lines( Result.unknown[,1],Result.unknown[,30], pch=8 , type="b", col="orange") # fEBVS
lines( Result.unknown[,1],Result.unknown[,34], pch=10, type="b", col="green" ) # Rebayes
lines( Result.unknown[,1],Result.unknown[,38], pch=19, type="b", col="red"   ) # FEBV
par(xpd=TRUE)
legend("topright",horiz=F, inset = c(-0.3,0.2), legend=c("Error Bar","Bonferroni","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","REBayes","F-EBV"),col=c("black","cyan3","blue","purple","green4","gold4", "brown","orange","green","red"),pch=c(NA,3,15,7,17,16,9,8,10,19), ncol=1, lty=1:1, cex=1)
dev.off()

## Length Ratio
figurename2 <- paste(path1,"LengRatio_",settname,".pdf",sep="")
pdf(figurename2,width=8,height=5)
par(mar=c(5,5,1,10))
Llim=max(Result.unknown[,11]/Result.unknown[,3],Result.unknown[,15]/Result.unknown[,3],Result.unknown[,35]/Result.unknown[,3])+0.1
if(Matrix==1)
  plot( Result.unknown[,1],Result.unknown[,3],xlab=TeX(r'($\tau^2/(1+\tau^2)$ )'),ylab="Length Ratio", ylim=c(0,Llim+0.2), xlim=c(0,1),col="white",main=bquote(paste("AR    ", ~rho,"=",.(rho),",   df=",.(df))), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
if(Matrix==2)
  plot( Result.unknown[,1],Result.unknown[,3],xlab=TeX(r'($\tau^2/(1+\tau^2)$ )'),ylab="Length Ratio", ylim=c(0,1),xlim=c(0,1),col="white",main=bquote(paste("Banded    ", ~rho,"=",.(rho),",   k=",.(k),",   df=",.(df))), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
if(Matrix==3)
  plot( Result.unknown[,1],Result.unknown[,3],xlab=TeX(r'($\tau^2/(1+\tau^2)$ )'),ylab="Length Ratio", ylim=c(0,Llim+0.2), xlim=c(0,1),col="white",main=bquote(paste("Sparse    ", ~rho," in (",.(a),",",.(b),"),   df=",.(df),",   Sparsity=",.(Sparsity),"%")), cex.lab=1.2, cex.axis=1.2, cex.main=1.2)


lines( Result.unknown[,1],Result.unknown[,11]/Result.unknown[,3], pch=15, type="b", col="blue"  ) # ELJS
lines( Result.unknown[,1],Result.unknown[,15]/Result.unknown[,3], pch=7 , type="b", col="purple") # TW
lines( Result.unknown[,1],Result.unknown[,19]/Result.unknown[,3], pch=17, type="b", col="green4") # Smyth
lines( Result.unknown[,1],Result.unknown[,23]/Result.unknown[,3], pch=16, type="b", col="gold4" ) # Vash
lines( Result.unknown[,1],Result.unknown[,27]/Result.unknown[,3], pch=9 , type="b", col="brown" ) # fEBV
lines( Result.unknown[,1],Result.unknown[,31]/Result.unknown[,3], pch=8 , type="b", col="orange") # fEBVS
lines( Result.unknown[,1],Result.unknown[,35]/Result.unknown[,3], pch=10, type="b", col="green" ) # ReBayes
lines( Result.unknown[,1],Result.unknown[,39]/Result.unknown[,3], pch=19, type="b", col="red"   ) # FEBV
##  legend("topright",horiz=F, inset = c(xloc,0), legend=c("ELJS","TW","Smyth","Vash","f-EBV","f-EBVS","F-EBV"),col=c("blue","purple","green4","gold4", "brown","orange","red"),pch=c(15,7,17,16,9,8,19), ncol=1, lty=1:1, cex=1)
dev.off()


  
