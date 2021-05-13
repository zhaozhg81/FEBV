xloc=-0.3

#-- Saving Figures ----------------------------------------------------------------------------------------------------------------------------

figurename1 <- paste(path1,"Coverage_",settname,".pdf",sep="")
pdf(figurename1,width=8,height=5)
par(mar=c(5,5,1,10) )
Llim=min(Result[,2],Result[,10],Result[,14],Result[,18],Result[,22])-0.1
 plot( Result[,1],Result[,2],xlab=TeX(r'($\tau^2$ )'),ylab="Coverage Prob", 
       ylim=c(Llim-0.1,1.05), xlim=c(0,3),col="white",
       main=bquote(paste("Leukemia Data ")), cex.lab=1.2, cex.axis=1.2, 
       cex.main=1.2)
lines( Result[,1],rep(0.95-qnorm(0.975)*sqrt(0.05*0.95/500),N), type="l", col="black",lwd=3,lty=2)
## lines( Result[,1],rep(0.95,N), type="l", col="black",lwd=3,lty=2)
lines( Result[,1],Result[,2] , pch=3 , type="b", col="cyan3" ) # Bonferroni
lines( Result[,1],Result[,10], pch=15, type="b", col="blue"  ) # ELJS
lines( Result[,1],Result[,14], pch=7 , type="b", col="purple") # TW
lines( Result[,1],Result[,18], pch=17, type="b", col="green4") # Smyth
lines( Result[,1],Result[,22], pch=16, type="b", col="gold4" ) # Vash
lines( Result[,1],Result[,26], pch=9 , type="b", col="brown" ) # fEBV
lines( Result[,1],Result[,30], pch=8 , type="b", col="orange") # fEBVS
lines( Result[,1],Result[,34], pch=10, type="b", col="green"   ) # rebayes
lines( Result[,1],Result[,38], pch=19, type="b", col="red"   ) # FEBV
par(xpd=TRUE)
legend("topright",horiz=F,inset=c(-0.3,0.2), legend=c("Error bar","Bonferroni","ELJS","TW","Smyth","Vash","f-EBV","f-EBVS",                                                      "REBayes","F-EBV"),
       col=c("black","cyan3","blue","purple","green4", "gold4","brown","orange","green","red"),
       pch=c(NA,3,15,7,17,16,9,8,10,19), ncol=1, lty=1:1, cex=1)

dev.off()

# Length Ratio
figurename2 <- paste(path1,"LengRatio_",settname,".pdf",sep="")
pdf(figurename2,width=8,height=5)
par(mar=c(5,5,1,10))
Llim=max(Result[,11]/Result[,3],Result[,15]/Result[,3],Result[,39]/Result[,3])
 plot( Result[,1],Result[,3],xlab=TeX(r'($\tau^2$ )'),ylab="Length Ratio", 
       ylim=c(0.1,Llim+0.01), xlim=c(0,3),col="white",
       main=bquote(paste("Leukemia Data ")), cex.lab=1.2, cex.axis=1.2, 
       cex.main=1.2)
lines( Result[,1],Result[,11]/Result[,3], pch=15, type="b", col="blue"  ) # ELJS
lines( Result[,1],Result[,15]/Result[,3], pch=7 , type="b", col="purple") # TW
lines( Result[,1],Result[,19]/Result[,3], pch=17, type="b", col="green4") # Smyth
lines( Result[,1],Result[,23]/Result[,3], pch=16, type="b", col="gold4" ) # Vash
lines( Result[,1],Result[,27]/Result[,3], pch=9 , type="b", col="brown" ) # fEBV
lines( Result[,1],Result[,31]/Result[,3], pch=8 , type="b", col="orange") # fEBVS
lines( Result[,1],Result[,35]/Result[,3], pch=10, type="b", col="green"   ) # rebayes
lines( Result[,1],Result[,39]/Result[,3], pch=19, type="b", col="red"   ) # FEBV

dev.off()

