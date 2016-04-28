#function to calculate abundance estimate and CV from dht with optional p_0 from separate surveys
#required because only some surveys were interobserver
#if p_0 and p_0_cv not specified function will return same output as dht
#don't use if there are coefficients in the detection function
#author: Lisa-Marie Harrison
#date: 27/04/2016

calcAbundanceAndCV<- function(det_fun, line_length, n_surveys, p_0=1, p_0_cv=0, group=TRUE) {
  
  #det_fun: detection function object from ddf
  #line_length: survey line length in km
  #n_surveys: number of surveys (needed or function deletes surveys with no observations)
  #p_0: p(0) estimate for first observer from separate interobserver trial surveys. Leave as default if no p_0 value available
  #p_0_cv: CV of p_0. Leave as default if no p_0_cv value available
  #group: boolean specifying if group CV or individual CV is required
  #returns: list containing abundance and CV
  
  dat <- det_fun$data
  p <- summary(det_fun)$average.p * p_0
  p_cv <- sqrt((summary(det_fun)$average.p.se / summary(det_fun)$average.p)^2 + p_0_cv^2)
  survey_area <- line_length *  det_fun$meta.data$width/1000
  
  if (group) {
    
    full_obs <- c(table(dat$Trial), rep(0, n_surveys - length(unique(dat$Trial))))
    var_n <- varn(rep(line_length, n_surveys), full_obs, type = "R2")
    n_L <- mean(full_obs/line_length)
    
    abund <- mean(full_obs)/(survey_area*p)*survey_area
    cv <- sqrt((sqrt(var_n)/n_L)^2 + p_cv^2)
    
  } else {
    
    n_inds <- c(aggregate(dat$size, by = list(dat$Trial), FUN = "sum")[, 2], rep(0, n_surveys - length(unique(dat$Trial))))
    var_n <- varn(rep(line_length, n_surveys), n_inds, type = "R2")
    n_L <- mean(n_inds/line_length)
    
    abund <- mean(n_inds)/(survey_area*p)*survey_area
    cv <- sqrt((sqrt(var_n)/n_L)^2 + p_cv^2)   
    
  }
  
  return (list(abund = abund, CV = cv))
  
}