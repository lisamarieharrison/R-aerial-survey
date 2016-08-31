#run multiple RJ-MCMC chains in runRjmcmc function in parallel and write output to separate files
#author: Lisa-Marie Harrison
#date: 25/07/2016

library(doParallel)
library(foreach)

source("C:/Users/43439535/Documents/Lisa/phd/aerial survey/R/R-aerial-survey/rjmcmc_parallel.R")

#choose number of parallel chains
n_chains <- 2 

#set up parallel back end
nodes <- detectCores() - 1
cl <- makeCluster(nodes)
registerDoParallel(cl)

#the ddf analysis found scale = 318 and shape = 3.6
#to convert ddf parameters exponentiate both and then add 1 to shape
foreach(chain = 1:n_chains, scale0 = c(100, 300), shape0 = c(1, 6), int0 = c(1, 2), species = rep("B", n_chains)) %dopar% {
  runRjmcmc(chain, scale0, shape0, int0, species)
}

stopCluster(cl)

