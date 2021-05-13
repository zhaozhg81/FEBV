This repo contains the functions to calculate the F-modeling based Empirical Bayes estimator of variances (F-EBV).
The function is stored in the file:
F_EBV.R

Usage:
FEB.est = FEB( sSq, df)

sSq: the sample variances
df: the degrees of freedom



The repo also contains all the code for running the simulation studies reported in the manuscript:

Yeil Kwon, Zhao, Z. (2019) On F-modelling based Empiricial Bayes Estimation of Variances. 

The aucillary files are stored in the directory R/
The source files for conducting the simulation studies are stored in the direcotry simulation/
The source files for the realdata analysis are stored in the directory realdata/

(i) Comparison of variance estimations post the selection
simulation/VarComp_[1-8].R

We have considered 4 different settings to generate the variances and there are two parameter
sets for each setting. When running the the simulation the first time, set PLOT as FALSE and
the simulation results will be stored under the folder Data/ and the figures are stored under the folder Figure/
If setting PLOT as TRUE, then the function will only plot the figures by loading the simulation results
stored in the directory Data/

(ii) Finite Bayes inference problem.
simulation/FinBay_[1-8].R

(iii) Comparison of the estimation of the mean parameters
simulation/MeanComp_[1-8].R

(iv) Comparison of the coverage probabilities of confidence intervals
simulation/CI_Ind_[1-8].R

(v) Concordance rate in the realdata analysis
realdata/run_dataanalysis.R

Warnings: when running the method REBayes, the function GVMix() might terminate and returns an error message.
This error is generated due to the optimization. If this happens, please rerun the code again.
