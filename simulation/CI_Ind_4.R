source("../R/F_Functions.R")
source("../R/F_Indep_Core.R")


#-- Setting Parameters -------------------------------------------------------------------------------------------------------------

prior   <- 2      # Prior:  1=IG, 2=Lognormal, 3=MixIG, 4=MixLognormal
a       <- -2
b       <- 1

##-- Setting Parameters ----------------------------------------------------------------------------------------------------------------------------------------------
      
G       <- 1000   # Dimension
numSim  <- 500    # The number of simulation
N       <- 40     # The Number of ticks for X-axis
n       <- 6      # Sample size  (df=n-1)

CI_Ind(G, numSim, N, n, prior, a, b, PLOT=FALSE, bw= 'ucv' )
