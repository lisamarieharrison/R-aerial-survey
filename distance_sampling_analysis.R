#analysis of aerial survey data using only Lisa's observations
#author: Lisa-Marie Harrison
#date: 22/10/2015\

setwd("C:/Users/Lisa/Documents/phd/aerial survey/R/data")
lisa_obs <- read.csv("lisa_full_observations.csv", header = T)
library(mrds)
library(chron)
library(knitr)

table(lisa_obs$Species)

lisa_obs$Trial <- as.numeric(lisa_obs$Date)

#total sightings
table(lisa_obs$Flight.Direction, lisa_obs$Species)

#combine all shark species into a single category, except hammerheads
lisa_obs$Species[lisa_obs$Species %in% c("BS", "W", "Wh", "S")] <- "S"


#---------------------------- DETECTION FUNCTIONS -----------------------------#



createData <- function(species, lisa_obs, truncate=NULL, direction) {
  
  lisa_obs <- lisa_obs[lisa_obs$Flight.Direction == direction & lisa_obs$Type == "S" 
                       & lisa_obs$Dist..from.transect <= truncate &
                       lisa_obs$Dist..from.transect !=0 & !is.na(lisa_obs$Dist..from.transect)
                       & lisa_obs$Species == species, ]

  total_observations <- cbind(1:nrow(lisa_obs), rep(1, nrow(lisa_obs)), rep(1, nrow(lisa_obs)), lisa_obs$Dist..from.transect, lisa_obs$Species, lisa_obs$Trial, lisa_obs$Number)
  total_observations <- data.frame(total_observations)
  colnames(total_observations) <- c("object", "observer", "detected", "distance", "species", "Trial", "size")

  if (species != "BOT") {
    total_observations <- total_observations[, 1:6]
  } else {
    total_observations <- total_observations[!is.na(total_observations$size), ]
  }
  
  return(total_observations)
  
}

total_observations <- createData(species = "S", lisa_obs, truncate = 1000, direction = "S")

p_total <- ddf(method = 'ds',dsmodel =~ cds(key = "gamma", formula=~1), 
               data = total_observations, meta.data = list(left = 50, width = 1000))
summary(p_total)
ddf.gof(p_total, main="Total observations goodness of fit")
plot(p_total, main = "Baitfish - NORTH")



area <- 265*0.3
obs.table <- cbind(rep(1, nrow(total_observations)), total_observations$Trial, total_observations)
colnames(obs.table)[1:2] <- c("Region.Label", "Sample.Label")

region.table <- data.frame(matrix(c(1, area), ncol = 2, byrow = T))
colnames(region.table) <- c("Region.Label", "Area")

sample.table <- data.frame(cbind(rep(1, length(unique(total_observations$Trial))), unique(total_observations$Trial), rep(265, length(unique(total_observations$Trial)))))
colnames(sample.table) <- c("Region.Label", "Sample.Label", "Effort")

d <- dht(p_total, region.table = region.table, sample.table = sample.table, obs.table = obs.table, 
         options = list(convert.units = 0.001))
d





#-------------------------- Environmental variables ---------------------------#

#remove single "LT" or "RT" or times will be incorrect
test <- dat.south$Type[dat.south$Type %in% c("LT", "RT")]
for (i in 2:length(test)) {
  if ((test[i] == "LT" & test[i - 1] == "LT") | (test[i] == "RT" & test[i - 1] == "RT")) {
    dat.south <- dat.south[-i, ]
  }
}  

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


dat.south <- lisa_obs[lisa_obs$Flight.Direction == "S", ]
dat.south <- dat.south[!dat.south$Time == "", ]

#sighting effort for southbound data
for (i in unique(dat.south$Date)) {
  for (j in 1:(nrow(dat.south[dat.south$Date == i, ]) - 1)) {
    
    if (dat.south[dat.south$Date == i, ]$Type[j] != "LT") {
      
      
      time.start <- chron(times. = as.character(dat.south$Time[dat.south$Date == i][j]), format = "h:m:s")
      time.stop  <- chron(times. = as.character(dat.south$Time[dat.south$Date == i][j + 1]), format = "h:m:s")
      mins <- hours(time.stop - time.start)*60 + minutes(time.stop - time.start)
      
      if (mins > 30) print(c(i, mins))
      
      #add minute differences to each environmental level
      for (k in c(22, 25, 26, 28)) {
        w <- which(colnames(envt.var.south) == as.name(paste(names(dat.south)[k], ".south.", dat.south[dat.south$Date == i, ][j, k], sep = "")))
        envt.var.south[1, w] <- envt.var.south[1, w] + mins
      }
    }
  }
}


#sea state
ss_hrs <- c(envt.var.south[1], envt.var.south[4], envt.var.south[5], envt.var.south[6])/60
t <- table(lisa_obs$Beaufort.Sea.State[lisa_obs$Species == "B" & lisa_obs$Flight.Direction == "S"])
sightings_hr <- matrix(c(t[1], t[2], t[3], t[4])/ss_hrs, ncol = 4) #wrong if no sightings at a level
colnames(sightings_hr) <- c("1", "2", "3", "4")
kable(sightings_hr, format = "pandoc", caption = "Number of sightings per hour effort at each sea state")


#turbidity
t_hrs <- c(envt.var.south[17:19])/60
t <- table(lisa_obs$Water.clarity[lisa_obs$Species == "S" & lisa_obs$Flight.Direction == "S"])
sightings_hr <- matrix(c(t[1], t[2], t[3])/t_hrs, ncol = 3) #wrong if no sightings at a level
colnames(sightings_hr) <- c("1", "2", "3")
kable(sightings_hr, format = "pandoc", caption = "Number of sightings per hour effort at each turbidity")


#sightings per survey at each wind strength  

wind_tab_B <- table(lisa_obs$Date[lisa_obs$Flight.Direction == "S" & lisa_obs$Species == "B"], 
                    lisa_obs$Wind.speed[lisa_obs$Flight.Direction == "S" & lisa_obs$Species == "B"])
wind_tab_B[wind_tab_B == 0] <- NA

wind_tab_BOT <- table(lisa_obs$Date[lisa_obs$Flight.Direction == "S" & lisa_obs$Species == "BOT"], 
                    lisa_obs$Wind.speed[lisa_obs$Flight.Direction == "S" & lisa_obs$Species == "BOT"])
wind_tab_BOT[wind_tab_BOT == 0] <- NA

wind_tab_S <- table(lisa_obs$Date[lisa_obs$Flight.Direction == "S" & lisa_obs$Species == "S"], 
                      lisa_obs$Wind.speed[lisa_obs$Flight.Direction == "S" & lisa_obs$Species == "S"])
wind_tab_S[wind_tab_S == 0] <- NA

plot(as.numeric(colnames(wind_tab_B)), colMeans(wind_tab_B, na.rm = TRUE), xlab = "Wind Speed (km/h)",
     ylab = "Mean number of sightings", type = "l", lwd = 2)
lines(as.numeric(colnames(wind_tab_BOT)), colMeans(wind_tab_BOT, na.rm = TRUE), col = "red", lwd = 2)
lines(as.numeric(colnames(wind_tab_S)), colMeans(wind_tab_S, na.rm = TRUE), col = "blue", lwd = 2)

legend("topright", c("Baitfish", "Dolphins", "Sharks"), col = c("black", "red", "blue"), lwd = 2, bty = "n")




