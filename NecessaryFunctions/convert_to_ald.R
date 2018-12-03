## Script to convert data to asymmetric Laplace evaluated at the maximum likelihood estimates.
#To use the data in MATLAB

##Author: Jami Jackson Mulgrave

library(ald)

setwd(getwd())

data <- matrix(unlist(read.csv(file="y_true_d.csv", header=FALSE, sep=",")), ncol = 1)

ald_MLE = mleALD(data)$par #Find the MLE of asymmetric Laplace using the parameters.

ald_cdf = pALD(data, mu = ald_MLE[1], sigma = ald_MLE[2], p = ald_MLE[3])

write.csv(ald_cdf, "ald_CDF_for_xmatrix.csv") #write the data for the cdf

write.csv(ald_MLE, "ald_MLE.csv") #write the MLE for ALD distribution