
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


#-- Showing Figures ----------------------------------------------------------------------------------------------------------------------------

par(mfrow=c(2,2))
barplot(log(L0[1,],10), ylab="log(MSE)", main=bquote(paste("[L0, Max]  ",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(L0[2,],10), ylab="log(MSE)", main=bquote(paste("[L0, 1%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(L0[3,],10), ylab="log(MSE)", main=bquote(paste("[L0, 5%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(L0[4,],10), ylab="log(MSE)", main=bquote(paste("[L0, All]  ",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)

barplot(log(L1[1,],10), ylab="log(MSE)", main=bquote(paste("[L1, Max]  ",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(L1[2,],10), ylab="log(MSE)", main=bquote(paste("[L1, 1%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(L1[3,],10), ylab="log(MSE)", main=bquote(paste("[L1, 5%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(L1[4,],10), ylab="log(MSE)", main=bquote(paste("[L1, All]  ",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)

barplot(log(L2[1,],10), ylab="log(MSE)", main=bquote(paste("[L2, Max]  ",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(L2[2,],10), ylab="log(MSE)", main=bquote(paste("[L2, 1%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(L2[3,],10), ylab="log(MSE)", main=bquote(paste("[L2, 5%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(L2[4,],10), ylab="log(MSE)", main=bquote(paste("[L2, All]  ",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)

barplot(log(Lp[1,],10), ylab="log(MSE)", main=bquote(paste("[L1p, Max]  ",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(Lp[2,],10), ylab="log(MSE)", main=bquote(paste("[L1p, 1%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(Lp[3,],10), ylab="log(MSE)", main=bquote(paste("[L1p, 5%]  " ,.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
barplot(log(Lp[4,],10), ylab="log(MSE)", main=bquote(paste("[L1p, All]  ",.(PRIOR), "    a=",.(a),",  b=",.(b))),name=namelist)
