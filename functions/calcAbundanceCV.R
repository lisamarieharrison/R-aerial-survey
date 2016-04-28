#function to calculate CV around abundance estimates from dht (just replicates dht CV calculator)
#author: Lisa-Marie Harrison
#date: 27/04/2016

calcAbundanceCV <- function(det_fun, line_length, n_surveys, group=TRUE) {
  
  #det_fun: detection function object from ddf
  #line_length: survey line length in km
  #n_surveys: number of surveys (needed or function deletes surveys with no observations)
  #group: boolean specifying if group CV or individual CV is required
  #returns: numeric CV
  
  dat <- det_fun$data
  p_cv <- summary(det_fun)$average.p.se / summary(det_fun)$average.p
  
  if (group) {
    
    full_obs <- c(table(dat$Trial), rep(0, n_surveys - length(unique(dat$Trial))))
    var_n <- varn(rep(line_length, n_surveys), full_obs, type = "R2")
    n_L <- mean(full_obs/line_length)
    
    cv <- sqrt((sqrt(var_n)/n_L)^2 + p_cv^2)
    
  } else {
    
    n_inds <- c(aggregate(dat$size, by = list(dat$Trial), FUN = "sum")[, 2], rep(0, n_surveys - length(unique(dat$Trial))))
    var_n <- varn(rep(line_length, n_surveys), n_inds, type = "R2")
    n_L <- mean(n_inds/line_length)
    
    cv <- sqrt((sqrt(var_n)/n_L)^2 + p_cv^2)   
    
  }
  
  return (cv)
  
}