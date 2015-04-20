#seasonal differences for aerial surveys
#author: Lisa-Marie Harrison
#date: 10/04/2015

setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/R/data")
source("C:/Users/Lisa/Documents/phd/aerial survey/R/code/R-aerial-survey/fun_calculate_envt_effort.R")
dat <- read.csv("full_data_set_20150410.csv", header = T)
library(chron)
library(stats)
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
dat_season_north$Species <- factor(dat_season_north$Species)
dat_season_south$Species <- factor(dat_season_south$Species)

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


#number of baitfish within 10km of each sighting
bot10k <- distance[bot, bait]
bot10k[bot10k > 10] <- NA
bot10k[bot10k <= 10] <- 1
fish_within_10k <- apply(bot10k, 1, sum, na.rm = TRUE)
hist(fish_within_10k, breaks = seq(0, max(fish_within_10k), by = 0.5))

white10k <- distance[white, bait]
white10k[white10k > 10] <- NA
white10k[white10k <= 10] <- 1
fish_within_10k <- apply(white10k, 1, sum, na.rm = TRUE)
hist(fish_within_10k, breaks = seq(0, max(fish_within_10k), by = 0.5))


#distance to baitfish from each sighting (not on same day) to check if sightings in high density prey areas
distance_all <- matrix(0, nrow = nrow(dat_season_north), ncol = nrow(dat_season_north))
for (i in 1:nrow(dat_season_north)) {
  for (j in 1:nrow(dat_season_north)) {
    distance_all[i, j] <- gcd.hf(deg2rad(dat_season_north$Lat[i]), deg2rad(dat_season_north$Long[i]), deg2rad(dat_season_north$Lat[j]), deg2rad(dat_season_north$Long[j]))    
  }
}

hist(distance_all[bot, bait])
hist(distance_all[white, bait])


#number of each species at each latitude
lat_tab <- table(dat$Species, round(dat$Lat, 2))
species <- which(rownames(lat_tab) %in% c("B", "BOT", "HH", "W"))

barplot(lat_tab[species[1], ], beside = TRUE, main = "# of baitfish sightings at each latitude", xlab = "latitude", ylab = "total sightings")
barplot(lat_tab[species[2], ], beside = TRUE, main = "# of dolphin sightings at each latitude", xlab = "latitude", ylab = "total sightings")
barplot(lat_tab[species[3], ], beside = TRUE, main = "# of hammerhead sightings at each latitude", xlab = "latitude", ylab = "total sightings")
barplot(lat_tab[species[4], ], beside = TRUE, main = "# of white shark sightings at each latitude", xlab = "latitude", ylab = "total sightings")


#tide data
tide_tab <- table(dat$Species, dat$Tide)
species <- which(rownames(lat_tab) %in% c("B", "BOT", "HH", "W"))
date_tab <- table(dat$Date, dat$Tide)
date_tab[date_tab > 0] <- 1
date_weight <- colSums(date_tab)

barplot(tide_tab[species[1], ]/date_weight, beside = TRUE, xlab = "tide height (m)", main = "baitfish sightings at tide heights")
barplot(tide_tab[species[2], ]/date_weight, beside = TRUE, xlab = "tide height (m)", main = "dolphin sightings at tide heights")
barplot(tide_tab[species[3], ]/date_weight, beside = TRUE, xlab = "tide height (m)", main = "hammerhead sightings at tide heights")
barplot(tide_tab[species[4], ]/date_weight, beside = TRUE, xlab = "tide height (m)", main = "white shark sightings at tide heights")


#swim direction
#no difference in N, S or C for flight directions
table(dat$Flight.Direction[dat$Species == "BOT"], dat$Swim.Direction[dat$Species == "BOT"])
table(dat$Flight.Direction[dat$Species == "W"], dat$Swim.Direction[dat$Species == "W"])
table(dat$Flight.Direction[dat$Species == "HH"], dat$Swim.Direction[dat$Species == "HH"])


#distance from transect
hist(dat_season_south$Dist..from.transect[dat_season_north$Species == "BOT"], xlab = "distance from transect (m)", main = "Southbound distance from transect - Bottlenose dolphin")

bait <- dat_season_south[dat_season_south$Species == "B" & dat_season_south$Length.Size %in% c("S", "M", "L"), ]
bait$Length.Size <- factor(bait$Length.Size)
boxplot(bait$Dist..from.transect ~ bait$Length.Size, xlab = "Baitfish size", ylab = "Distance from transect", main = "Southbound distance from transect - baitfish")


# -------------- ARE SIGHTINGS CLUSTERED ALONG THE TRANSECT? ---------------- #

#Hopkin's statistic
calculateHopkins <- function(species, dat, direction) {
  
  #function to calculate Hopkin's statistic for clustering
  #formula from Skaug (2006)
  #species = string denoting species
  #dat = full data set
  #direction = string denoting flight direction (N or S)
  
  lat <- seq(min(na.omit(dat$Lat)), max(na.omit(dat$Lat)), length.out = 10000)
  

  species_obs <- dat[dat$Species == species & !is.na(dat$Lat) & dat$Flight.Direction == direction, ]
  n <- nrow(species_obs)
  
  
  chosen_point <- sample(lat, n, replace = FALSE)
  
  Xi <- 0
  for(i in 1:length(chosen_point)) {
    w <- which.min(abs(species_obs$Lat - chosen_point[i]))
    Xi[i] <- gcd.hf(deg2rad(chosen_point[i]), 0, deg2rad(species_obs$Lat[w]), 0)    
  }
  
  Yi <- 0
  for(i in 1:n) {
    rem_obs <- species_obs$Lat[-i]
    w <- which.min(abs(rem_obs - species_obs$Lat[i]))
    Yi[i] <- gcd.hf(deg2rad(species_obs$Lat[i]), 0, deg2rad(rem_obs[w]), 0)    
  }
  
  Hf <- sum(Xi^2)/sum(Yi^2) #Hf > 0.5 indicates clustering
  
  return(Hf)
  
}

calculateHopkins(species = "BOT", dat = dat, direction = "N")


#kolmogorov-smirnoff to compare poisson
calculateKolmogorovSmirnoff <- function(species, dat, direction) {
 
  #tests if sighting rate is poisson using Kolmogorov-Smirnoff test
  #species = string containing species code
  #dat = full data set
  #direction = string containing flight direction
  
  lat <- seq(min(na.omit(dat$Lat)), max(na.omit(dat$Lat)), length.out = 10000)
  species_obs <- dat[dat$Species == species & !is.na(dat$Lat) & dat$Flight.Direction == direction, ]
  n <- nrow(species_obs)  
  chosen_point <- sample(lat, n, replace = FALSE)
  
  event_spacing <- 0
  for (i in 2:length(species_obs$Lat)) {
    event_spacing[i] <- gcd.hf(deg2rad(sort(species_obs$Lat)[i]), 0, deg2rad(sort(species_obs$Lat)[i - 1]), 0)
  }
  
  sim_event_spacing <- 0
  for (i in 2:length(chosen_point)) {
    sim_event_spacing[i] <- gcd.hf(deg2rad(sort(chosen_point)[i]), 0, deg2rad(sort(chosen_point)[i - 1]), 0)
  }
  
  lambda = mean(sim_event_spacing)
  pois_rgen <- rpois(n, lambda)
  
  ks.test(event_spacing, pois_rgen)    
  
}

calculateKolmogorovSmirnoff("B", dat, "S")


# ----------------------- TEST SEASONAL DIFFERENCES -------------------------- #
season_tab <- table(dat_season_north$Date, dat_season_north$Season)
season_tab[season_tab > 0] <- 1
days_per_season <- colSums(season_tab)
table(dat_season_north$Species, dat_season_north$Season)/rep(days_per_season, each = 4)


d <- dat_season_north
d$Species <- rep(1, nrow(d))
sighting_tab <- table(d$Date, d$Species)

lm.dat <- as.matrix(cbind(rownames(sighting_tab), sighting_tab[, 1], as.character(dat$Season[dat$Type == "SS" & dat$Flight.Direction == "N"])))
colnames(lm.dat) <- c("date", "sightings", "season")

lm.season <- lm(as.numeric(lm.dat[, 2]) ~ lm.dat[, 3])
summary(lm.season)

boxplot(as.numeric(lm.dat[, 2]) ~ lm.dat[, 3])





