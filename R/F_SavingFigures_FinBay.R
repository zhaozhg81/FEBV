
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

par(mfrow=c(2,2))
namelist=c("sSq",	"ELJS",	"TW",	"Smyth",	"Vash",	"fEBV",	"fEBVS",	"FEBV")

figurename1 <- paste(path1,"FinBay_",settname,"_L0.pdf",sep="")
pdf(figurename1,width=8,height=5)
barplot(log(L0,10), ylab="log(MSE)", main=bquote(paste("Finite Bayes:",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
dev.off()

figurename1 <- paste(path1,"FinBay_",settname,"_L1.pdf",sep="")
pdf(figurename1,width=8,height=5)
barplot(log(L1,10), ylab="log(MSE)", main=bquote(paste("Finite Byaes:" ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
dev.off()

figurename1 <- paste(path1,"FinBay_",settname,"_L2.pdf",sep="")
pdf(figurename1,width=8,height=5)
barplot(log(L2,10), ylab="log(MSE)", main=bquote(paste("Finite Bayes:" ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
dev.off()

figurename1 <- paste(path1,"FinBay_",settname,"_Lp.pdf",sep="")
pdf(figurename1,width=8,height=5)
barplot(log(Lp,10), ylab="log(MSE)", main=bquote(paste("Finite Bayes:",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
dev.off()


#-- Showing Figures ----------------------------------------------------------------------------------------------------------------------------

par(mfrow=c(2,2))
barplot(log(L0,10), ylab="MSE (log_scale)", main=bquote(paste("[L0, Max]  ",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(L1,10), ylab="MSE (log_scale)", main=bquote(paste("[L1, 1%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(L2,10), ylab="MSE (log_scale)", main=bquote(paste("[L2, 5%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(Lp,10), ylab="MSE (log_scale)", main=bquote(paste("[Lp, All]  ",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
