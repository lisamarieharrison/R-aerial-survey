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

#the ddf analysis found scale = 318 and shape = 3.6 for dolphins, and 197 and 303 for fish
#to convert ddf parameters exponentiate both and then add 1 to shape
foreach(chain = 1:n_chains, scale0 = c(180, 230), shape0 = c(4, 3), int0 = c(0.2, 1), species = rep("B", n_chains), truncate_left = c(0, 0), truncate_right = c(1000, 1000)) %dopar% {
  runRjmcmc(chain, scale0, shape0, int0, species, truncate_left, truncate_right)
}

stopCluster(cl)

