#find effort at each environmental variable
#Author: Lisa-Marie Harrison
#Date: 07/11/2014

setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/R/data")
dat <- read.csv("environmental_effort.csv", header = T)[1:2427, ]
library(chron)

dat$Time <- chron(times. = dat$Time, format = "h:m:s")

dat.north <- dat[dat$Flight.Direction == "N", ]
dat.south <- dat[dat$Flight.Direction == "S", ]

#------------------------ EFFORT FOR NORTH BOUND FLIGHTS ----------------------#

#initialize environmental variables
envt.var.north <- matrix(0, ncol = 27)
colnames(envt.var.north) <- c(
      "Beaufort.Sea.State.north.1",
      "Beaufort.Sea.State.north.0",
      "Beaufort.Sea.State.north.1.5",
      "Beaufort.Sea.State.north.2",
      "Beaufort.Sea.State.north.2.5",
      "Beaufort.Sea.State.north.3",
      "Beaufort.Sea.State.north.4", 
      "Cloud.cover.north.0",
      "Cloud.cover.north.1",
      "Cloud.cover.north.2",
      "Cloud.cover.north.3",
      "Cloud.cover.north.4",
      "Cloud.cover.north.5",
      "Cloud.cover.north.6",
      "Cloud.cover.north.7",
      "Cloud.cover.north.8",
      "Water.clarity.north.1",
      "Water.clarity.north.2",
      "Water.clarity.north.3",
      "Glare.north.0",
      "Glare.north.10",
      "Glare.north.30",
      "Glare.north.15",
      "Glare.north.25",
      "Glare.north.5",
      "Glare.north.20",
      "Glare.north.1"
      )


#sighting effort for northbound data
for (i in unique(dat.north$Date)) {
  for (j in 1:(nrow(dat.north[dat.north$Date == i, ]) - 1)) {
    
    if (dat.north[dat.north$Date == i, ]$Type[j] != "LT") {
      
      time.start <- dat.north$Time[dat.north$Date == i][j] 
      time.stop  <- dat.north$Time[dat.north$Date == i][j + 1]
      mins <- hours(time.stop - time.start)*60 + minutes(time.stop - time.start)
      
      #add minute differences to each environmental level
      for (k in c(5, 8, 9, 11)) {
        w <- which(colnames(envt.var.north) == as.name(paste(names(dat.north)[k], ".north.", dat.north[dat.north$Date == i, ][j, k], sep = "")))
        envt.var.north[1, w] <- envt.var.north[1, w] + mins
      }
    }
  }
}



#------------------------ EFFORT FOR SOUTH BOUND FLIGHTS ----------------------#

#initialize environmental variables
envt.var.south <- matrix(0, ncol = 30)
colnames(envt.var.south) <- c(
  "Beaufort.Sea.State.south.1",
  "Beaufort.Sea.State.south.0",  
  "Beaufort.Sea.State.south.1.5",
  "Beaufort.Sea.State.south.2",
  "Beaufort.Sea.State.south.3",
  "Beaufort.Sea.State.south.4", 
  "Beaufort.Sea.State.south.5", 
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
for (i in unique(dat.south$Date)) {
  for (j in 1:(nrow(dat.south[dat.south$Date == i, ]) - 1)) {
    
    if (dat.south[dat.south$Date == i, ]$Type[j] != "LT") {
      
      
      time.start <- dat.south$Time[dat.south$Date == i][j] 
      time.stop  <- dat.south$Time[dat.south$Date == i][j + 1]
      mins <- hours(time.stop - time.start)*60 + minutes(time.stop - time.start)
      
      
      #add minute differences to each environmental level
      for (k in c(5, 8, 9, 11)) {
        w <- which(colnames(envt.var.south) == as.name(paste(names(dat.south)[k], ".south.", dat.south[dat.south$Date == i, ][j, k], sep = "")))
        envt.var.south[1, w] <- envt.var.south[1, w] + mins
      }
    }
  }
}

envt.mins <- cbind(envt.var.north, envt.var.south)
colnames(envt.mins) <- c(colnames(envt.var.north), colnames(envt.var.south))

write.csv(envt.mins, "environmental_effort_minutes.csv", row.names = F)



