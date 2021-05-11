PLOT=TRUE
library(REBayes)
library(latex2exp)

    
options(digits=3)

##-- Functions ----------------------------------------------------------------------------------------

source("../../R/F_Functions.R")

##-- Create Correlation Matrix ------------------------------------------------------------------------

G       <- 1000   # Dimension
numSim  <- 500   # The number of simulation
N       <- 20     # The Number of ticks for X-axis
df      <- n-1    # Degree of Freedom
    Matrix  <- 3

idx        <- sample.int(G*G, size =G*G*prop, replace = FALSE)
elem       <- a+b*runif(G*G*prop)
elemM      <- array(0,c(G*G,1))
elemM[idx] <- elem
A          <- matrix(elemM,G,G)
A1         <- A-diag(diag(A))+diag(G)
SM         <- A1%*%t(A1)
B          <- diag(1/sqrt((diag(SM))))
R          <- B%*%SM%*%B
B          <- t(chol(R))
Sparsity   <- length(c(R[abs(R)==0]))/G/G*100
quant      <- quantile(c(R[abs(R)>0 & abs(R)<1]),c(seq(0,1,0.05)))

if(PLOT==FALSE)
  {
    ##-- Common Core --------------------------------------------------------------------------------------
    source("../../R/F_CommonCore.R")
    
    ##-- Saving RData -------------------------------------------------------------------------------------------------------------------
    
    path      <- paste("./Data/")
    path1     <- paste("./Figure/")
    settname  <- paste("Sparse_a_",a,"_b_",b,"_df_",df,"_prop_",prop,sep="")
    filename  <- paste(path, settname, ".Rdata", sep="")
    save(Result.unknown, file = filename)
  }else{
    ##-- Loading RData -------------------------------------------------------------------------------------------------------------------
    
    path      <- paste("./Data/")
    path1     <- paste("./Figure/")
    settname  <- paste("Sparse_a_",a,"_b_",b,"_df_",df,"_prop_",prop,sep="")
    filename  <- paste(path, settname, ".Rdata", sep="")
    load(filename)
  }

##-- Saving Figures -----------------------------------------------------------------------------------


path1     <- paste("./Figure/")
source("../../R/F_SavingFigures_Depend.R")

path1     <-"/home/zhaozhg/Dropbox/Apps/Overleaf/On F-modelling based Empirical Bayes Estimation of Variances/figure/"
source("../../R/F_SavingFigures_Depend.R")

##-----------------------------------------------------------------------------------------------------
