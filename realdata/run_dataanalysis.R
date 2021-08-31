source("../realdata/Concordance_Colon.R")
source("../realdata/Concordance_Leukemia.R")

## concordance_colon(PLOT=TRUE)
## concordance_leukemia(PLOT=TRUE)

## To rerun the code, do the following
concordance_colon(numSim=500, PLOT=FALSE, bw='defaults')
concordance_leukemia(numSim=500, PLOT=FALSE, bw='defaults')


