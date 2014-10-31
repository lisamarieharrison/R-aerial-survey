#aerial survey exploratory analysis
#includes data from October 2013, Summer 13/14 and April 2014
#27/08/2014

dat <- read.csv(file = "C:/Users/Lisa/Documents/phd/aerial survey/data/sightings.csv", header = T)


#number of sightings of each species
table(dat$Species)

#histograms by major species of distance from transect for Northbound surveys
par(mfrow = c(2, 2))
hist(dat$Dist..from.transect[dat$Flight.Direction == "N" & dat$Species == "B"], 
     xlab = "distance from transect", main = "Baitfish Schools")
hist(dat$Dist..from.transect[dat$Flight.Direction == "N" & dat$Species == "BOT"], 
     xlab = "distance from transect", main = "Bottlenose Dolphin Pods")
hist(dat$Dist..from.transect[dat$Flight.Direction == "N" & dat$Species == "HH"], 
     xlab = "distance from transect", main = "Hammerhead Sharks")
hist(dat$Dist..from.transect[dat$Flight.Direction == "N" & dat$Species == "W"], 
     xlab = "distance from transect", main = "White Sharks")


#histograms by major species of distance from transect for Southbound surveys
par(mfrow = c(2, 2))
hist(dat$Dist..from.transect[dat$Flight.Direction == "S" & dat$Species == "B"], 
     xlab = "distance from transect", main = "Baitfish Schools")
hist(dat$Dist..from.transect[dat$Flight.Direction == "S" & dat$Species == "BOT"], 
     xlab = "distance from transect", main = "Bottlenose Dolphin Pods")
hist(dat$Dist..from.transect[dat$Flight.Direction == "S" & dat$Species == "HH"], 
     xlab = "distance from transect", main = "Hammerhead Sharks")
hist(dat$Dist..from.transect[dat$Flight.Direction == "S" & dat$Species == "W"], 
     xlab = "distance from transect", main = "White Sharks")






