calcEnvtEffort <- function(dat) {

  #function to calculate the number of hours at each environmental variable level
  #dat = matrix of all waypoints
  #return = matrix of minutes at each environmental variable level
  
  library(chron)
  
  #initialize environmental variables
  envt.var.south <- matrix(0, ncol = 27)
  colnames(envt.var.south) <- c(
    "Beaufort.Sea.State.south.1",
    "Beaufort.Sea.State.south.2",
    "Beaufort.Sea.State.south.3",
    "Beaufort.Sea.State.south.4", 
    "Cloud.cover.south.0",
    "Cloud.cover.south.1",
    "Cloud.cover.south.2",
    "Cloud.cover.south.3",
    "Cloud.cover.south.4",
    "Cloud.cover.south.5",
    "Cloud.cover.south.6",
    "Cloud.cover.south.7",
    "Cloud.cover.south.8",
    "Water.clarity.south.1",
    "Water.clarity.south.2",
    "Water.clarity.south.3",
    "Glare.south.0",
    "Glare.south.10",
    "Glare.south.30",
    "Glare.south.15",
    "Glare.south.25",
    "Glare.south.5",
    "Glare.south.20",
    "Glare.south.35",
    "Glare.south.50",
    "Glare.south.60",
    "Glare.south.100"
  )
  
  #sighting effort for southbound data
  for (i in unique(dat$Date)) {
    for (j in 1:(nrow(dat[dat$Date == i, ]) - 1)) {
      
      if (dat[dat$Date == i, ]$Type[j] != "LT") {
        
        
        time.start <- chron(times. = as.character(dat$Time[dat$Date == i][j]), format = "h:m:s")
        time.stop  <- chron(times. = as.character(dat$Time[dat$Date == i][j + 1]), format = "h:m:s")
        mins       <- hours(time.stop - time.start)*60 + minutes(time.stop - time.start)
        
        if (mins > 30) {
          message(paste(i, mins)) #check for times that have been typed in incorrectly
        }
        
        #add minute differences to each environmental level
        for (k in c(24, 27, 28, 30)) {
          w <- which(colnames(envt.var.south) == as.name(paste(names(dat)[k], ".south.", dat[dat$Date == i, ][j, k], sep = "")))
          envt.var.south[1, w] <- envt.var.south[1, w] + mins
        }
      }
    }
  }
  
  return(envt.var.south)
  
}