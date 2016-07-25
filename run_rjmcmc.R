library(doParallel)
library(foreach)

source("C:/Users/43439535/Documents/Lisa/phd/aerial survey/R/R-aerial-survey/rjmcmc_parallel.R")

#choose number of parallel chains
n_chains <- 2 

#set up parallel back end
nodes <- detectCores() - 1
cl <- makeCluster(nodes)
registerDoParallel(cl)

foreach(chain=1:n_chains, scale0=c(100, 200), shape0=c(1, 10), int0 = c(1, 5)) %dopar% {
  runRjmcmc(chain, scale0, shape0, int0)
}

stopCluster(cl)



