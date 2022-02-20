source("../R/F_Functions.R")
source("../R/F_Finite_Bayes.R")
G <- 1000
numSim <- 500
n <- 6
a <- 3
b <- 1
prior <- 6
PLOT=FALSE

F_FinBayes(G, numSim, n, prior, a, b, PLOT)


## #-- Setting Parameters -------------------------------------------------------------------------------------------------------------

## prior   <- 1      # Prior:  1=IG, 2=Lognormal, 3=Gamma, 4=IGMix, 5=GamMix
## a       <- 3
## b       <- 1



