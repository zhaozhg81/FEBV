CI_AR <- function(G, numSim, N, n, prior, a, b, rho, PLOT=FALSE, bw='ucv')
{

  
  options(digits=3)
  
 
  ##-- Create Correlation Matrix --------------------------------------------------------------------------
  
  ## G       <- 1000   # Dimension
  ## numSim  <- 500   # The number of simulation
  ## N       <- 20     # The Number of ticks for X-axis

  df      <- n-1    # Degree of Freedom
  Matrix  <- 1
  
  R=diag(G)
  for(ii in 1:(G-1))
    {
      for(jj in (ii+1):(G))
        {
          R[ii,jj]=(rho)^(jj-1)
        }
    }
  R=R+t(R)-diag(G)
  B=t(chol(R))
  
  ##-- Common Core ----------------------------------------------------------------------------------------
  
  if(PLOT==FALSE)
    {
      source("../R/F_CommonCore.R")
      
      ##-- Saving RData -------------------------------------------------------------------------------------------------------------------
      
      path      <- paste("./Data/")
      path1     <- paste("./Figure/")
      settname  <- paste("AR_rho_",rho,"_df_",df,sep="")
      filename  <- paste(path, settname, ".Rdata", sep="")
      save(Result.unknown, file = filename)
    }else{
      path      <- paste("./Data/")
      path1     <- paste("./Figure/")
      settname  <- paste("AR_rho_",rho,"_df_",df,sep="")
      filename  <- paste(path, settname, ".Rdata", sep="")
      load(filename)    
    }
  
  ##-- Saving Figures ------------------------------------------------------------------------------------
  path1     <- paste("./Figure/")
  source("../F_SavingFigures_Depend.R")
  
  ## path1     <-"/home/zhaozhg/Dropbox/Apps/Overleaf/On F-modelling based Empirical Bayes Estimation of Variances/figure/"
  ## source("../F_SavingFigures_Depend.R")
  
  ##------------------------------------------------------------------------------------------------------
}



