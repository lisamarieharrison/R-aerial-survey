setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/R/data")
data <- read.csv(file = "season aerial survey.csv", header = T, fill = T)
dist <- read.csv(file = "season distance offshore.csv", header = F, fill = T)
envt <- read.csv(file = "sighting envt.csv", header = T, fill = T)
oct  <- read.csv(file = "oct data.csv", header = T, fill = T)
attach(data)

names(data)
table(Flight.Direction, Species)
table(Species)

dat <- subset(data, Species %in% c("B", "BOT", "HH", "W"))

#total number of individuals
sum(na.omit(Number[Flight.Direction == "S" & Species == "W"]))

#number of white sharks at Swansea-Blacksmiths beach
sum(Number[Location == "Swansea-Blacksmiths" & Species == "W"])

bb <- subset(dat, Species == "B")
bot <- subset(dat, Species == "BOT")
hh <- subset(dat, Species == "HH")
white <- subset(dat, Species == "W")

par(mfrow = c(1, 4))
hist(bb$Dist..from.transect, main = "baitfish")
hist(bot$Dist..from.transect, main = "bottlenose dolphin")
hist(hh$Dist..from.transect, main = "hammerhead")
hist(white$Dist..from.transect, main = "white shark")

#distance offshore by species



boxplot(dat$Dist..from.transect ~ dat$Species, xlab = "Species", ylab = "Distance offshore (m)", xaxt = 'n')
axis(side = 1, at = 1:8, c("", "Baitfish", "Bottlenose dolphin", "", "Hammerhead", "", "", "White Shark"))
title("Distance offshore by species for Southbound transects")


#sightings by cloud cover
a <- table(Cloud.cover, Species)

plot(a[, 2], type = "l", ylim = c(0, 45), xlab = "cloud cover", ylab = "number of sightings", lwd = 2.5)
points(a[, 3], type = "l", col = "red", lwd = 2.5)
points(a[, 4], type = "l", col = "blue", lwd = 2.5)
points(a[, 9], type = "l", col = "orange", lwd = 2.5)
title("Number of sightings by cloud cover")
legend("topleft", c("Baitfish", "Bottlenose Dolphin", "Hammerhead Shark", "White Shark"), col = c("black", "red", "blue", "orange")
       , lwd = 2.5, bty = "n")



#sightings by sea state
a <- table(Beaufort.Sea.State, Species)

plot(a[, 2], type = "l", ylim = c(0, 130), xlab = "Sea state", ylab = "number of sightings", xaxt = "n", lwd = 2.5)
axis(side = 1, at = 1:4, c(1, 2, 3, 4))
points(a[, 3], type = "l", col = "red", lwd = 2.5)
points(a[, 4], type = "l", col = "blue", lwd = 2.5)
points(a[, 9], type = "l", col = "orange", lwd = 2.5)
title("Number of sightings by sea state")
legend("topright", c("Baitfish", "Bottlenose Dolphin", "Hammerhead Shark", "White Shark"), col = c("black", "red", "blue", "orange")
       , lwd = 2.5, bty = "n")





#sightings by turbidity
a <- table(Water.clarity, Species)

plot(a[, 2], type = "l", ylim = c(0, 180), xlab = "Water clarity", ylab = "number of sightings", xaxt = "n", lwd = 2.5)
axis(side = 1, at = 1:3, c("Excellent", "Good", "Poor"))
points(a[, 3], type = "l", col = "red", lwd = 2.5)
points(a[, 4], type = "l", col = "blue", lwd = 2.5)
points(a[, 9], type = "l", col = "orange", lwd = 2.5)
title("Number of sightings by water clarity")
legend("topright", c("Baitfish", "Bottlenose Dolphin", "Hammerhead Shark", "White Shark"), col = c("black", "red", "blue", "orange")
       , lwd = 2.5, bty = "n")


#sighting rate for each species
a <- table(Date[Flight.Direction == "N"], Species[Flight.Direction == "N"])
b <- table(Date[Flight.Direction == "S"], Species[Flight.Direction == "S"])

boxplot(a[, 1], b[, 1], a[, 2], b[, 2], a[, 3], b[, 3], a[, 4], b[, 4], xaxt = "n", xlab = "Species", ylab = "Sightings per transect"
        , border = c("black", "red"), lwd = 2.5)
axis(side = 1, at = c(1.5, 3.5, 5.5, 7.5), c("Baitfish", "Bottlenose dolphin", "Humpback Whale", "White Shark"))
title("Sighting rate for north and southbound transects")
legend("topright", c("Northbound", "Southbound"), col = c("black", "red"), lwd = 2.5, bty = "n")



