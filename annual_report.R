#analysis for 2013/2014 aerial survey annual report
#author: Lisa-Marie Harrison
#date: 05/11/2014

setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/R/data")
dat <- read.csv("annual_data.csv", header = T)
av.ind <- read.csv("annual_average_individuals.csv", header = T)
int <- read.csv("annual_interobserver.csv", header = T)[1:491, ]

wind.dir <- rep("OFF", nrow(dat))
wind.dir[dat$Wind.direction >= 0 & dat$Wind.direction <= 180] <- "ON"
wind.dir[is.na(dat$Wind.direction)] <- NA

#plot average number of sightings per day across seasons by species
bar <- table(dat$season, dat$Species)
bar[1, ] <- round(bar[1, ]/6)
bar[2, ] <- round(bar[2, ]/8)
bar[3, ] <- round(bar[3, ]/13)
par(mar=c(10,4,2,2))
barplot(bar, beside = T, legend = c("autumn (n = 6)", "spring (n = 8)", "summer (n = 13)"), names.arg = c("baitfish", "bottlenose dolphin", "large individual fish", "humpback whale", "hammerhead shark", "pinniped", "ray", "unidentified shark", "turtle", "white shark"), 
        cex.names = 1.2, ylab = "Mean sightings/day", las = 2, args.legend = list(x = 48, y = 45, bty = "n"))

#plot average number of individuals per sighting
barplot(as.matrix(av.ind[, 2:ncol(av.ind)]), beside = T, legend = c("autumn", "spring", "summer"), 
        cex.names = 1.2, las = 2, args.legend = list(x = 45, y = 14, bty = "n"), names.arg = c("baitfish", "bottlenose dolphin", "large individual fish", "humpback whale", "hammerhead shark", "pinniped", "ray", "unidentified shark", "turtle", "white shark"), 
        ylab = "Average individuals/sighting")


#-------------------------- ENVIRONMENTAL CONDITIONS --------------------------#

#plot of distance from transect by flight direction for each sea state
par(mar = c(8, 4, 2, 2))
at <- c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14)
lvl1 <- c("N", "S")
lvl2 <- c("1 (n = 125)", "1.5 (n = 33)", "2 (n = 814)", "3 (n = 119)", "4 (n = 14)")
boxplot(dat$Dist..from.transect ~ dat$Flight.Direction + dat$Beaufort.Sea.State, at = at, xaxt = "n", xlab = "", las = 2, ylab = "Sighting distance from transect (m)")
axis(1, at = at, labels = rep(lvl1, 5), tick = FALSE)
mtext(lvl2, 1, line = 3, at = c(1.5, 4.5, 7.5, 10.5, 13.5), font = 2)
mtext("Beaufort Sea State (bold) and Flight Direction (North/South)", 1, line = 5)

#plot of distance from transect by flight direction by wind direction (offshore/onshore)
par(mar = c(8, 4, 2, 2))
at <- c(1:4)
lvl1 <- c("Offshore", "Onshore")
lvl2 <- c("North", "South")
boxplot(dat$Dist..from.transect ~ dat$Flight.Direction + wind.dir, at = at, xaxt = "n", xlab = "", las = 2, ylab = "Sighting distance from transect (m)")
axis(1, at = at, labels = rep(lvl1, 2), tick = FALSE)
mtext(lvl2, 1, line = 3, at = c(1.5, 3.5), font = 2)
mtext("Wind direction", 1, line = 5)


#plot of distance from transect by flight direction by wind speed (kmph)
plot(dat$Wind.speed, dat$Dist..from.transect, xlab = "Wind speed (km/h)", las = 2, ylab = "Sighting distance from transect (m)", pch = 19)

#plot of distance from transect by flight direction by cloud cover
par(mar = c(8, 4, 2, 2))
at <- c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 19, 20, 22, 23, 25, 26)
lvl1 <- c("N", "S")
lvl2 <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")
boxplot(dat$Dist..from.transect ~ dat$Cloud.cover + dat$Flight.Direction, at = at, xaxt = "n", xlab = "", las = 2, ylab = "Sighting distance from transect (m)")
axis(1, at = at, labels = rep(lvl1, 9), tick = FALSE)
mtext(lvl2, 1, line = 3, at = seq(1.5, 26.5, by = 3), font = 2)
mtext("Cloud Cover (0-8)", 1, line = 5)


#lm to test distance from transect in different weather conditions
dat.north <- dat[dat$Flight.Direction == "N" & is.na(dat$Dist..from.transect) == FALSE, ]
dat.north$wind.dir <-  wind.dir[dat$Flight.Direction == "N" & is.na(dat$Dist..from.transect) == FALSE]
north.lm <- lm(log(dat.north$Dist..from.transect) ~ dat.north$Beaufort.Sea.State + dat.north$wind.dir + dat.north$Wind.speed + dat.north$Cloud.cover + dat.north$Water.clarity + dat.north$Glare)
summary(north.lm)
 
dat.south <- dat[dat$Flight.Direction == "S" & is.na(dat$Dist..from.transect) == FALSE, ]
dat.south$wind.dir <-  wind.dir[dat$Flight.Direction == "S" & is.na(dat$Dist..from.transect) == FALSE]
south.lm <- lm(log(dat.south$Dist..from.transect) ~ dat.south$Beaufort.Sea.State + dat.south$wind.dir + dat.south$Wind.speed + dat.south$Cloud.cover + dat.south$Water.clarity + dat.south$Glare)
summary(south.lm)


#plot of number of sightings per hour in each environmental condition level
par(mfrow = c(2, 2))

sightings_ss_hour <- table(dat$Flight.Direction, dat$Beaufort.Sea.State)
sightings_ss_hour[1, ] <- sightings_ss_hour[1, ]/c(191/60, 74/60, 1431/60, 311/60, 1)
sightings_ss_hour[2, ] <- sightings_ss_hour[2, ]/c(193/60, 23/60, 1403/60, 444/60, 58/60)

barplot(sightings_ss_hour, beside = T, xlab = "Beaufort Sea State", ylab = "sightings/hour effort")


sightings_cc_hour <- table(dat$Flight.Direction, dat$Cloud.cover)
sightings_cc_hour[1, ] <- sightings_cc_hour[1, ]/c(499, 394, 85, 141, 191, 182, 54, 96, 372)*60
sightings_cc_hour[2, ] <- sightings_cc_hour[2, ]/c(425, 289, 131, 415, 326, 138, 8, 222, 169)*60

barplot(sightings_cc_hour, beside = T, xlab = "Cloud cover", ylab = "sightings/hour effort")

sightings_wc_hour <- table(dat$Flight.Direction, dat$Water.clarity)
sightings_wc_hour[1, ] <- sightings_wc_hour[1, ]/c(120, 1462, 432)*60
sightings_wc_hour[2, ] <- sightings_wc_hour[2, ]/c(76, 1696, 351)*60

barplot(sightings_wc_hour, beside = T, xlab = "Water clarity (1 = excellent, 2 = good, 3 = poor)", ylab = "sightings/hour effort")

sightings_g_hour <- table(dat$Flight.Direction, dat$Glare)
sightings_g_hour[1, ] <- sightings_g_hour[1, ]/c(1117, 283, 508, 74, 1, 1, 11, 1, 1)*60
sightings_g_hour[2, ] <- sightings_g_hour[2, ]/c(1462, 161, 182, 106, 48, 43, 4, 45, 38)*60 

barplot(sightings_g_hour, beside = T, xlab = "Glare (%)", ylab = "sightings/hour effort")


#inter-observer sightings
table(int$Observer, int$Species)






