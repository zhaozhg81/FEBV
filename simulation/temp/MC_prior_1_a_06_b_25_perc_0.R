source("../R/F_Functions.R")
source("../R/F_mean_performance.R")


#-- Setting Parameters -------------------------------------------------------------------------------------------------------------

prior   <- 1       # Prior:  1=IG, 2=Lognormal, 3=mix IG, 4= mix lognormal
a       <- 6
b       <- 2.5
perc    <- 0.001   # Proportion of the number of selected parameters
N <- 20

G <- 1000
numSim <- 100
n <- 6
PLOT=FALSE


F_mean_Comparison(G, numSim, N, n, prior, a, b, perc, bw='defaults', PLOT=FALSE)

