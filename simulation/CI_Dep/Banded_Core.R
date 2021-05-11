PLOT=TRUE
library(REBayes)
library(latex2exp)


options(digits=3)

##-- Functions ------------------------------------------------------------------------------------------

source("../../R/F_Functions.R")

##-- Create Correlation Matrix --------------------------------------------------------------------------

G       <- 1000   # Dimension
numSim  <- 500  # The number of simulation
N       <- 20     # The Number of ticks for X-axis
df      <- n-1    # Degree of Freedom
Matrix  <- 2


R=diag(G);
for(ii in 1:(G-1))
  {
    for(jj in seq((ii+1),min(G,(ii+k)),1))
      {
        R[ii,jj]=rho
      }
  }
R=R+t(R)-diag(G)
B=t(chol(R))

##-- Common Core ----------------------------------------------------------------------------------------

if(PLOT==FALSE)
  {
    source("../../R/F_CommonCore.R")
    
    ##-- Saving RData -------------------------------------------------------------------------------------------------------------------
    
    path      <- paste("./Data/")
    path1     <- paste("./Figure/")
    settname  <- paste("Banded_rho_",rho,"_k_",k,"_df_",df,sep="")
    filename  <- paste(path, settname, ".Rdata", sep="")
    save(Result.unknown, file = filename)
  }else{
   
    ##-- Loading RData -------------------------------------------------------------------------------------------------------------------
    
    path      <- paste("./Data/")
    path1     <- paste("./Figure/")
    settname  <- paste("Banded_rho_",rho,"_k_",k,"_df_",df,sep="")
    filename  <- paste(path, settname, ".Rdata", sep="")
    load(filename)
  }        
    
    ##-- Saving Figures -------------------------------------------------------------------------------------

    path1     <- paste("./Figure/")
    source("../../R/F_SavingFigures_Depend.R")
    
    path1     <-"/home/zhaozhg/Dropbox/Apps/Overleaf/On F-modelling based Empirical Bayes Estimation of Variances/figure/"
    source("../../R/F_SavingFigures_Depend.R")
    
    ##-------------------------------------------------------------------------------------------------------

