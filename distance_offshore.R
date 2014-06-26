setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/data")
data <- read.table(file = "distance_offshore.txt", header = T, fill = T)
attach(data)

table(Species)

plot(Species, Dist, xlab = "Species: B = Baitfish (n = 34), BOT = Bottlenose (n = 30), HB = Humpback (n = 2), S = Shark (n = 1), W = White Shark (n = 10)"
     , ylab = "Distance offshore (m)")
title("Number of sightings for each species by distance offshore (Southbound transects)")


par(mfrow = c(2, 2))
hist(Dist[Species == "B"], xlab = "Distance offshore (m)", main = "Sightings by distance offshore for Baitfish (n = 34)", breaks = seq(1, 600, 25))
hist(Dist[Species == "BOT"], xlab = "Distance offshore (m)", main = "Sightings by distance offshore for Bottlenose dolphins (n = 30)", breaks = seq(1, 600, 25))
hist(Dist[Species == "HB"], xlab = "Distance offshore (m)", main = "Sightings by distance offshore for Humpback Whales (n = 2)", breaks = seq(1, 600, 25))
hist(Dist[Species == "W"], xlab = "Distance offshore (m)", main = "Sightings by distance offshore for White Sharks (n = 10)", breaks = seq(1, 600, 25))

data <- read.table(file = "distance_offshore_n.txt", header = T, fill = T)
attach(data)

table(Species)

par(mfrow = c(1, 1))
plot(Species, Dist, xlab = "Species: B = Baitfish (n = 21), BOT = Bottlenose (n = 13), HB = Humpback (n = 2), W = White Shark (n = 5)"
     , ylab = "Distance from transect (m)")
title("Number of sightings for each species by distance from transect (Northbound transects)")


par(mfrow = c(2, 2))
hist(Dist[Species == "B"], xlab = "Distance from transect (m)", main = "Sightings by distance from transect for Baitfish (n = 21)", breaks = seq(1, 525, 25))
hist(Dist[Species == "BOT"], xlab = "Distance from transect (m)", main = "Sightings by distance from transect for Bottlenose dolphins (n = 13)", breaks = seq(1, 525, 25))
hist(Dist[Species == "HB"], xlab = "Distance from transect (m)", main = "Sightings by distance from transecte for Humpback Whales (n = 2)", breaks = seq(1, 525, 25))
hist(Dist[Species == "W"], xlab = "Distance from transect (m)", main = "Sightings by distance from transect for White Sharks (n = 5)", breaks = seq(1, 550, 25))



data <- read.table(file = "distance_from_transect.txt", header = T, fill = T)
attach(data)

table(Species)

par(mfrow = c(1, 1))
plot(Species, Dist, xlab = "Species: B = Baitfish (n = 55), BOT = Bottlenose (n = 43), HB = Humpback (n = 4), S = shark (n = 1), W = White Shark (n = 15)"
     , ylab = "Distance from transect (m)")
title("Number of sightings for each species by distance from transect")


par(mfrow = c(2, 2))
hist(Dist[Species == "B"], xlab = "Distance from transect (m)", main = "Sightings by distance from transect for Baitfish (n = 55)", breaks = seq(1, 600, 25))
hist(Dist[Species == "BOT"], xlab = "Distance from transect (m)", main = "Sightings by distance from transect for Bottlenose dolphins (n = 43)", breaks = seq(1, 600, 25))
hist(Dist[Species == "HB"], xlab = "Distance from transect (m)", main = "Sightings by distance from transecte for Humpback Whales (n = 4)", breaks = seq(1, 600, 25))
hist(Dist[Species == "W"], xlab = "Distance from transect (m)", main = "Sightings by distance from transect for White Sharks (n = 15)", breaks = seq(1, 600, 25))



