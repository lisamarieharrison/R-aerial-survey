library(plyr)
library(doParallel)

nodes <- detectCores() - 1
cl <- makeCluster(nodes)
registerDoParallel(cl)

l_ply(1:2, runRjmcmc, .parallel = TRUE)

stopCluster(cl)



