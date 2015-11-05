#inter-observer aerial survey  analysis for 22 Lisa-Vic southbound flights
#author: Lisa-Marie Harrison
#date: 20/10/2015


dat <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/interobserver_20151105.csv", header = T)
library(knitr)
library(mrds)

#source required functions
source("C:/Users/Lisa/Documents/phd/aerial survey/R/code/R-aerial-survey/functions/checkSameSighting.R")
source("C:/Users/Lisa/Documents/phd/aerial survey/R/code/R-aerial-survey/functions/findOverlappingObservations.R")
source("C:/Users/Lisa/Documents/phd/aerial survey/R/code/R-aerial-survey/functions/createDistanceDataInterobserver.R")
source("C:/Users/Lisa/Documents/phd/aerial survey/R/code/R-aerial-survey/functions/calcAbundanceInterobserver.R")


dat <- dat[dat$Type == "S", ]
dat <- dat[dat$Species != "", ]
dat$Species[dat$Species == "B "] <- "B"
dat$Species <- factor(dat$Species)


kable(table(dat$Observer, dat$Species), format = "pandoc", caption = "Total sightings by each observer")

  
lisa <- dat[dat$Observer == "Lisa", ]
vic  <- dat[dat$Observer == "Vic", ]

#separate observations into overlap (both seen) and missed observations
obs <- findOverlappingObservations(lisa, vic)
overlap     <- obs$overlap
lisa_missed <- obs$obs1_missed
vic_missed  <- obs$obs2_missed

#table of sightings
obs_numbers <- cbind(nrow(lisa_missed), nrow(vic_missed), nrow(overlap))
colnames(obs_numbers) <- c("Vic only", "Lisa only", "Both seen")
kable(obs_numbers, format = "pandoc", caption = "Number of sightings seen by each observer")

#missed sightings by species  
par(mfrow = c(1, 1))  
total_missed <- rbind(lisa_missed, vic_missed)
missed_table <- table(total_missed$Observer, total_missed$Species)
x <- barplot(missed_table, beside = TRUE, legend = c("Vic", "Lisa"), main = "Number of missed sightings by species", xlab = "Species")

#number of missed sightings per survey
par(mfrow = c(1, 2))
hist(table(total_missed$Observer, total_missed$Date)[1, ], main = "Vic missed", xlab = "", xlim = c(0, 35), col = "lightgrey")
hist(table(total_missed$Observer, total_missed$Date)[2, ], main = "Lisa missed", xlab = "", xlim = c(0, 35), col = "lightgrey")

#missed sightings by date
total_missed$Date <- chron(dates. = as.character(total_missed$Date), format = "d/m/y")
par(mfrow = c(1, 1))
plot(table(total_missed$Observer, total_missed$Date)[1, ], type = "l", ylim = c(0, 35), 
     xlab = "Date", xaxt = "n", ylab = "# of missed sightings", lwd = 2)
points(table(total_missed$Observer, total_missed$Date)[2, ], col = "red", type = "l", lwd = 2)
legend("topleft", c("Vic", "Lisa"), col = c("black", "red"), lwd = 2, bty = "n", cex = 0.8)
axis(side = 1, at = 1:22, colnames(table(total_missed$Observer, total_missed$Date)))
title("Missed sightings by date")

#total number of sightings by each observer

lisa_missed <- poolAllSharks(lisa_missed)
overlap     <- poolAllSharks(overlap)
vic_missed  <- poolAllSharks(vic_missed)

vic_seen  <- table(lisa_missed$Species) + table(overlap$Species)
lisa_seen <- table(vic_missed$Species) + table(overlap$Species)

#vic percentage of Lisa
percent_of_lisa <- vic_seen/lisa_seen

#---------------------------- DETECTION FUNCTIONS -----------------------------#

total_observations <- createDistanceDataInterobserver(1, overlap, lisa_missed, vic_missed, truncate = 1000)

p_total <- ddf(method="io", mrmodel =~ glm(~distance), dsmodel =~ cds(key = "gamma", formula=~1),
               data = total_observations, meta.data = list(left = 50, width = 1000, point = FALSE))


summary(p_total)
#ddf.gof(p_total, main="Total observations goodness of fit")
#plot(p_total)


d <- calcAbundanceInterobserver(265, dataframe = total_observations, p_total, type = p_total$method)
d




