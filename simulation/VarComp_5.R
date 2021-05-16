source("../R/F_Functions.R")
source("../R/F_variance_performance.R")
G <- 1000
numSim <- 500
n <- 6
a <- 10
b <- 1
prior <- 3
PLOT=TRUE
bw='ucv'

F_variance_Comparison(G, numSim, n, prior, a, b, bw, PLOT)


## #-- Setting Parameters -------------------------------------------------------------------------------------------------------------

## prior   <- 1      # Prior:  1=IG, 2=Lognormal, 3=Gamma, 4=IGMix, 5=GamMix
## a       <- 3
## b       <- 1



