#inter-observer aerial survey preliminary analysis
setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/R/data")
dat <- read.csv("interobserver_20150407.csv", header = T)
library(knitr)


#total sightings seen by each observer
table(dat$Observer, dat$Species)

#number of overlapping sightings
lisa <- dat[dat$Observer == "Lisa", ]
vic <- dat[dat$Observer == "Vic", ]


checkSame <- function(r1, r2) {
  #checks if 2 observations are the same and returns boolean
  if (r1$Date == r2$Date & substr(r1$Time, 1, 5) == substr(r2$Time, 1, 5) & r1$Species == r2$Species) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}



findOverlap <- function(set1, set2, i, j,overlap, lisa_missed, vic_missed) {
  
  
  if (j > nrow(set2)) {
    print(paste(i, "No overlapping observation found. Try next observation"))
    lisa_missed <- rbind(lisa_missed, set1[i, ])
    return(list(overlap = overlap, vic_missed = vic_missed, lisa_missed = lisa_missed))
  }
  
  if (checkSame(set1[i, ], set2[j, ])) {
    overlap <- rbind(overlap, set1[i, ])
    print("Success: found a same observation!")
    vic_missed <- vic_missed[-which(vic_missed$Date == set1$Date[i] & substr(vic_missed$Time, 1, 5) == substr(set1$Time[i], 1, 5) & vic_missed$Species == set1$Species[i])[1], ]
    return(list(overlap = overlap, vic_missed = vic_missed, lisa_missed = lisa_missed))
  } else {
    findOverlap(set1, set2, i = i, j = j + 1, overlap = overlap, lisa_missed = lisa_missed, 
                vic_missed = vic_missed)
  }
  
}

overlap <- dat[0, ]
vic_missed <- lisa
lisa_missed <- dat[0, ]

for (i in 1:nrow(vic)) {
  set2 <- vic_missed[vic_missed$Date == vic$Date[i], ]
  
  result <- findOverlap(set1 = vic, set2 = set2, i = i, j = 1, 
                        overlap = overlap, lisa_missed = lisa_missed, 
                        vic_missed = vic_missed)
  overlap <- result$overlap
  lisa_missed <- result$lisa_missed
  vic_missed <- result$vic_missed
}
lisa_missed <- lisa_missed[-1, ]


#plot missing sightings by observer and species

total_missed <- rbind(lisa_missed, vic_missed)
missed_table <- table(total_missed$Observer, total_missed$Species)
x = barplot(missed_table, beside = TRUE, legend = c("Vic", "Lisa"), main = "Number of missed sightings by species", xlab = "Species")

par(mfrow = c(1, 2))
hist(table(total_missed$Observer, total_missed$Date)[1, ], main = "Vic missed", xlab = "", xlim = c(0, 35), col = "lightgrey")
hist(table(total_missed$Observer, total_missed$Date)[2, ], main = "Lisa missed", xlab = "", xlim = c(0, 35), col = "lightgrey")


#------------------------ ENVIRONMENTAL EFFORT ANALYSIS --------------------------#

dat <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/environmental_effort_20150407.csv", header = T)
library(chron)

dat$Time <- chron(times. = dat$Time, format = "h:m:s")

dat.south <- dat[dat$Flight.Direction == "S", ]

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

ss_percent <- c(envt.var.south[1], envt.var.south[4], envt.var.south[5], envt.var.south[6])/60
t <- table(total_missed$Beaufort.Sea.State)
missed_hr <- matrix(c(t[1]+t[2], t[3], t[4], t[5])/ss_percent, ncol = 4)
colnames(missed_hr) <- c("1", "2", "3", "4")
kable(missed_hr, format = "pandoc", caption = "Number of missed sightings per hour effort at each sea state")

total_missed$Date <- chron(dates. = total_missed$Date, format = "d/m/y")

plot(table(total_missed$Observer, total_missed$Date)[1, ], type = "l", ylim = c(0, 35), 
     xlab = "Date", xaxt = "n", ylab = "# of missed sightings", lwd = 2)
points(table(total_missed$Observer, total_missed$Date)[2, ], col = "red", type = "l", lwd = 2)
legend("topleft", c("Vic", "Lisa"), col = c("black", "red"), lwd = 2, bty = "n", cex = 0.8)
axis(side = 1, at = 1:16, colnames(table(total_missed$Observer, total_missed$Date)))

