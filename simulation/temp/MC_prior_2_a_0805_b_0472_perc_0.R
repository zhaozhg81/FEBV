source("../R/F_Functions.R")
source("../R/F_mean_performance.R")

#-- Setting Parameters -------------------------------------------------------------------------------------------------------------

prior   <-  2      # Prior:  1=IG, 2=Lognormal, 3=Gamma, 4=IGMix, 5=GamMix
a       <- -0.805
b       <-  0.472
perc    <-  0.001    # Proportion of the number of selected parameters
N <- 20


G <- 1000
numSim <- 100
n <- 6
PLOT=FALSE


F_mean_Comparison(G, numSim, N, n, prior, a, b, PLOT=FALSE)
