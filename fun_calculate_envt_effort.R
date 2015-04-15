#function to calculate environmental effort
#author: Lisa-Marie Harrison
#date: 10/04/2015


calculateEffort <- function(dat) {
  
  library(chron)
  
  #initialize environmental variables
  envt.var <- matrix(0, ncol = 27)
  colnames(envt.var) <- c(
    "Beaufort.Sea.State.1",
    "Beaufort.Sea.State.0",
    "Beaufort.Sea.State.1.5",
    "Beaufort.Sea.State.2",
    "Beaufort.Sea.State.2.5",
    "Beaufort.Sea.State.3",
    "Beaufort.Sea.State.4", 
    "Cloud.cover.0",
    "Cloud.cover.1",
    "Cloud.cover.2",
    "Cloud.cover.3",
    "Cloud.cover.4",
    "Cloud.cover.5",
    "Cloud.cover.6",
    "Cloud.cover.7",
    "Cloud.cover.8",
    "Water.clarity.1",
    "Water.clarity.2",
    "Water.clarity.3",
    "Glare.0",
    "Glare.10",
    "Glare.30",
    "Glare.15",
    "Glare.25",
    "Glare.5",
    "Glare.20",
    "Glare.1"
  )
  
  
  #sighting effort for bound data
  for (i in unique(dat$Date)) {
    for (j in 1:(nrow(dat[dat$Date == i, ]) - 1)) {
      
      if (dat[dat$Date == i, ]$Type[j] != "LT") {
        
        time.start <- dat$Time[dat$Date == i][j] 
        time.stop  <- dat$Time[dat$Date == i][j + 1]
        mins <- hours(time.stop - time.start)*60 + minutes(time.stop - time.start)
        
        #add minute differences to each environmental level
        for (k in c(22, 25, 26, 29)) {
          w <- which(colnames(envt.var) == as.name(paste(names(dat)[k], ".", dat[dat$Date == i, ][j, k], sep = "")))
          envt.var[1, w] <- envt.var[1, w] + mins
        }
        if(mins > 30) print(c(i, mins))
      }
    }
  }
  
  total <- sum(envt.var[1:7])
  
  
  return(list(envt.var = envt.var, total = total))
  
}

