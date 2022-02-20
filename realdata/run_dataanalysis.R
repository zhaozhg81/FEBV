source("../realdata/Concordance_Colon.R")
source("../realdata/Concordance_Leukemia.R")
source("../realdata/Prostate.R")

## concordance_colon(PLOT=TRUE)
## concordance_leukemia(PLOT=TRUE)

## To rerun the code, do the following
concordance_colon(numSim=500, PLOT=TRUE, bw='defaults')
concordance_leukemia(numSim=500, PLOT=TRUE, bw='defaults')


## concordance_prostate(numSim=500, PLOT=FALSE, bw='defaults')
