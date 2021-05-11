source("../R/F_Functions.R")
source("../R/F_Indep_Core.R")


#-- Setting Parameters -------------------------------------------------------------------------------------------------------------

prior   <- 1      # Prior:  1=IG, 2=Lognormal, 3=Gamma, 4=IGMix, 5=GamMix
a       <- 6
b       <- 2.5

##-- Setting Parameters ----------------------------------------------------------------------------------------------------------------------------------------------
      
G       <- 1000   # Dimension
numSim  <- 500    # The number of simulation
N       <- 20     # The Number of ticks for X-axis
n       <- 6      # Sample size  (df=n-1)

CI_Ind(G, numSim, N, n, prior, a, b, PLOT=FALSE, bw= 'defaults' )
