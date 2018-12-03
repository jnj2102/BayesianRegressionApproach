########################
#Code to run the BDMCMC method using the BDGraph package
#For the Simulation section of my paper
#Bayesian Regression Approach
#Author: Jami Jackson Mulgrave
#
########################

current_dir <- getwd()
setwd(current_dir)  #set the working directory

#clear the workspace
rm(list = ls())

library(R.matlab)
library(BDgraph)
library(Matrix)


source("BDGraph_copula.R") #call the function


reps <- 100

###########Simulation combination: n=100, p=50, sparsity = AR2#######
result_list <- list()

data <- readMat('CholeskyDecomp_p50_n100_AR2_prior_FFN.mat', package = "R.matlab") #Read in all the data


for (iters in 1:reps) {
  cat('Iteration = ', iters)

set.seed(100 + iters)
  
#Run this in a loop for each iteration 
  
  #read in the Sigma_true 
  
  sigma_true <- data$sigma.true
  

  #read in the omega_true
  
  omega_true <- data$omega.true
  
  p = data$p
  
  #read in the x matrix (observed data)
  
xmat <- data$x.matrix.n100.p50[[iters]][[1]]

#nonparanormal truncation with graphical lasso
result_list[[iters]] <- BDGraph_copula(omega_true, sigma_true, xmat,p)


}

###################################################################

#Save the data for latex later

save.image(file = "BDGraph_p50_n100_AR2.rdata")
     
