#seasonal differences for aerial surveys
#author: Lisa-Marie Harrison
#date: 10/04/2015

setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/R/data")
source("C:/Users/Lisa/Documents/phd/aerial survey/R/code/R-aerial-survey/fun_calculate_envt_effort.R")
dat <- read.csv("full_data_set_20150410.csv", header = T)
library(chron)
dat$Time <- chron(times. = dat$Time, format = "h:m:s")
dat$Date <- chron(dates. = as.character(dat$Date), format = "d/m/y")

dat$Beaufort.Sea.State[dat$Beaufort.Sea.State == 1.5] <- 1
dat_north <- dat[dat$Flight.Direction == "N", ]
dat_south <- dat[dat$Flight.Direction == "S", ]

#calculate number of hours flown in each direction
effort_north <- calculateEffort(dat_north)
envt_effort <- effort_north$envt.var
hours_flown_north <- effort_north$total


effort_south <- calculateEffort(dat_south)
envt_effort_south <- effort_south$envt.var
hours_flown_south <- effort_south$total

# ------------------------- SPECIES HISTOGRAMS BY DATE ----------------------- #

dat <- dat[order(dat$Date), ]
sightings_tab <- table(dat$Date, dat$Species)
d <- range(dat$Date)

#baitfish
x <- seq(d[1], d[2])
y <- rep(0, length(x))
y[x %in% unique(dat$Date)] <- sightings_tab[, which(colnames(sightings_tab) == "B")]
barplot(y, main = "Baitfish sightings", col = "gray23")
axis(1, at = c(10, 150, 350, 450, 560), c("Summer 13/14", "April 14", "Sept 14", "Summer 14/15", "April 15"))

#bottlenose dolphin
x <- seq(d[1], d[2])
y <- rep(0, length(x))
y[x %in% unique(dat$Date)] <- sightings_tab[, which(colnames(sightings_tab) == "BOT")]
barplot(y, main = "Bottlenose dolphin sightings", col = "gray23")
axis(1, at = c(10, 150, 350, 450, 560), c("Summer 13/14", "April 14", "Sept 14", "Summer 14/15", "April 15"))

#hammerhead shark
x <- seq(d[1], d[2])
y <- rep(0, length(x))
y[x %in% unique(dat$Date)] <- sightings_tab[, which(colnames(sightings_tab) == "HH")]
barplot(y, main = "Hammerhead shark sightings", col = "gray23")
axis(1, at = c(10, 150, 350, 450, 560), c("Summer 13/14", "April 14", "Sept 14", "Summer 14/15", "April 15"))

#White shark 
x <- seq(d[1], d[2])
y <- rep(0, length(x))
y[x %in% unique(dat$Date)] <- sightings_tab[, which(colnames(sightings_tab) == "W")]
barplot(y, main = "White shark sightings", col = "gray23")
axis(1, at = c(10, 150, 350, 450, 560), c("Summer 13/14", "April 14", "Sept 14", "Summer 14/15", "April 15"))



# -------------------- SIGHTINGS AT EACH ENVIRONTMENTAL LEVEL ---------------- #
dat_season_north <- dat_north[dat_north$Species %in% c("B", "BOT", "HH", "W"), ]
dat_season_south <- dat_south[dat_south$Species %in% c("B", "BOT", "HH", "W"), ]


#sea state
north_tab <- table(dat_season_north$Beaufort.Sea.State, dat_season_north$Species)
species <- which(colnames(north_tab) %in% c("B", "BOT", "HH", "W"))
north_tab_effort <- north_tab[, species]/(envt_effort[c(1, 4, 6)]/60)

barplot(north_tab_effort, beside = T, legend = T, main = "Sightings/hr at each sea state - north", xaxt = "n")
axis(1, at = c(2, 6, 10, 14), c("Baitfish", "Bottlenose Dolphin", "Hammerhead Shark", "White shark"))

#cloud cover
north_tab <- table(dat_season_north$Cloud.cover, dat_season_north$Species)
north_tab_effort <- north_tab[, species]/(envt_effort[c(8:16)]/60)

barplot(north_tab_effort, beside = T, main = "Sightings/hr at each cloud cover - north", xaxt = "n", legend = T, args.legend = c(cex = 0.5))
axis(1, at = c(5, 18, 29, 40), c("Baitfish", "Bottlenose Dolphin", "Hammerhead Shark", "White shark"))


#turbidity
north_tab <- table(dat_season_north$Water.clarity, dat_season_north$Species)
north_tab_effort <- north_tab[, species]/(envt_effort[c(17:19)]/60)

barplot(north_tab_effort, beside = T, main = "Sightings/hr at each turbidity - north", xaxt = "n", legend = T)
axis(1, at = c(2, 6, 10, 14), c("Baitfish", "Bottlenose Dolphin", "Hammerhead Shark", "White shark"))



#sea state
south_tab <- table(dat_season_south$Beaufort.Sea.State, dat_season_south$Species)
species <- which(colnames(south_tab) %in% c("B", "BOT", "HH", "W"))
south_tab_effort <- south_tab[1:3, species]/(envt_effort_south[c(1, 4, 6)]/60)

barplot(south_tab_effort, beside = T, legend = T, main = "Sightings/hr at each sea state - south", xaxt = "n")
axis(1, at = c(2, 6, 10, 14), c("Baitfish", "Bottlenose Dolphin", "Hammerhead Shark", "White shark"))

#cloud cover
south_tab <- table(dat_season_south$Cloud.cover, dat_season_south$Species)
south_tab_effort <- south_tab[, species]/(envt_effort_south[c(8:16)]/60)

barplot(south_tab_effort, beside = T, main = "Sightings/hr at each cloud cover - south", xaxt = "n", legend = T, args.legend = c(cex = 0.5))
axis(1, at = c(5, 18, 29, 40), c("Baitfish", "Bottlenose Dolphin", "Hammerhead Shark", "White shark"))


#turbidity
south_tab <- table(dat_season_south$Water.clarity, dat_season_south$Species)
south_tab_effort <- south_tab[1:3, species]/(envt_effort_south[c(17:19)]/60)

barplot(south_tab_effort, beside = T, main = "Sightings/hr at each turbidity - south", xaxt = "n", legend = T)
axis(1, at = c(2, 6, 10, 14), c("Baitfish", "Bottlenose Dolphin", "Hammerhead Shark", "White shark"))


# ----------------------- DISTANCE BETWEEN SIGHTINGS ------------------------- #

deg2rad <- function(deg) {
  #converts degrees to radians
  #input: degree coordinate
  #returns: radian coordinate 
  
  return(deg*pi/180)
}

gcd.hf <- function(lat1, long1, lat2, long2) {
  #calculates distance between two coordinates using the Haversine formula (hf)
  #input: radian latitude and longitude coordinates
  #returns: distance between coordinates in km
  
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat  <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1, sqrt(a)))
  d = R * c
  return(d) 
  
}

#calculate distances between all sightings on the same day
distance <- matrix(0, nrow = nrow(dat_season_north), ncol = nrow(dat_season_north))
for (i in 1:nrow(dat_season_north)) {
  for (j in 1:nrow(dat_season_north)) {
    if (dat_season_north$Date[i] == dat_season_north$Date[j]) {
      distance[i, j] <- gcd.hf(deg2rad(dat_season_north$Lat[i]), deg2rad(dat_season_north$Long[i]), deg2rad(dat_season_north$Lat[j]), deg2rad(dat_season_north$Long[j]))    
    }
  }
}

species <- dat_season_north$Species
bot <- which(species == "BOT")
bait <- which(species == "B")
white <- which(species == "W")
distance[distance == 0] <- NA
hist(distance[bot, bait], xlab = "km", main = "Distance between bottlenose dolphin and baitfish sightings", col = "lightgrey")
hist(distance[white, bait], xlab = "km", main = "Distance between white sharks and baitfish sightings", col = "lightgrey")


#minimum distance to fish from each sighting
min_dist_to_fish <- apply(distance[bot, bait], 1, min, na.rm = TRUE)
min_dist_to_fish[min_dist_to_fish == "Inf"] <- NA
hist(min_dist_to_fish, breaks = 150, main = "Distance to closest baitfish - bottlenose dolphin")


min_dist_to_fish <- apply(distance[white, bait], 1, min, na.rm = TRUE)
min_dist_to_fish[min_dist_to_fish == "Inf"] <- NA
hist(min_dist_to_fish, breaks = 150, main = "Distance to closest baitfish - white shark")







