poolAllSharks <- function(dat) {
  
  #pools all species of sharks into a single categoty, "S", except hammerheads
  #dat = data.frame of observations
  
  #if S level doesn't already exist, add it
  if (!"S" %in% levels(dat$Species)) {
    levels(dat$Species) <- c(levels(dat$Species), "S")
  }
  
  levels(dat$Species) <- c(levels(dat$Species), "S")
  
  dat$Species[dat$Species %in% c("BS", "W", "Wh")] <- "S"
  dat$Species <- factor(dat$Species)
  
  return(dat)
  
}