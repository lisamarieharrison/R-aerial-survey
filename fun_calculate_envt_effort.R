#function to calculate environmental effort
#author: Lisa-Marie Harrison
#date: 10/04/2015


calculateEffort <- function(dat) {
  
  library(chron)
  
  #initialize environmental variables
  names <- c(
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
    "Glare.1", 
    paste(rep("Wind.speed", max(na.omit(dat[, 25])) - min(na.omit(dat[, 25])) + 1), 
          rep(".", max(na.omit(dat[, 25])) - min(na.omit(dat[, 25])) + 1), 
          c(seq(min(na.omit(dat[, 25])), max(na.omit(dat[, 25])))), sep = "")
  )
  envt.var <- matrix(0, ncol = length(names))
  colnames(envt.var) <- names
  
  #sighting effort for bound data
  for (i in unique(dat$Date)) {
    for (j in 1:(nrow(dat[dat$Date == i, ]) - 1)) {
      
      if (dat[dat$Date == i, ]$Type[j] != "LT") {
        
        time.start <- dat$Time[dat$Date == i][j] 
        time.stop  <- dat$Time[dat$Date == i][j + 1]
        mins <- hours(time.stop - time.start)*60 + minutes(time.stop - time.start)
        
        #add minute differences to each environmental level
        for (k in c(23, 25, 26, 27, 29)) {
          w <- which(colnames(envt.var) == as.name(paste(names(dat)[k], ".", dat[dat$Date == i, ][j, k], sep = "")))
          envt.var[1, w] <- envt.var[1, w] + mins
        }
        if(mins > 30) print(c(as.character(unique(dat$Date)[unique(dat$Date) == i]), as.character(time.start), as.character(time.stop), mins))
      }
    }
  }
  
  total <- sum(envt.var[1:7])
  
  
  return(list(envt.var = envt.var, total = total))
  
}

