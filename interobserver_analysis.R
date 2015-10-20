#inter-observer aerial survey preliminary analysis
dat <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/interobserver_20151019.csv", header = T)
library(knitr)
library(mrds)

dat <- dat[dat$Type == "S", ]
dat <- dat[dat$Species != "", ]
dat$Species[dat$Species == "B "] <- "B"
dat$Species <- factor(dat$Species)

dat$Trial <- as.numeric(dat$Date)

kable(table(dat$Observer, dat$Species), format = "pandoc", caption = "Total sightings by each observer")

  
lisa <- dat[dat$Observer == "Lisa", ]
vic  <- dat[dat$Observer == "Vic", ]


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
    #print(paste(i, "No overlapping observation found. Try next observation"))
    lisa_missed <- rbind(lisa_missed, set1[i, ])
    return(list(overlap = overlap, vic_missed = vic_missed, lisa_missed = lisa_missed))
  }
  
  if (checkSame(set1[i, ], set2[j, ])) {
    overlap <- rbind(overlap, set1[i, ])
    #print("Success: found a same observation!")
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


#---------------------------- DETECTION FUNCTIONS -----------------------------#

#detection functions - hazard rate chosen as best using AIC

total_observations <- rbind(overlap, lisa_missed, vic_missed)
total_observations$Species[total_observations$Species %in% c("BS", "W", "Wh")] <- "S"
total_observations <- total_observations[total_observations$Species == "B", ]
total_observations <- total_observations[total_observations$Dist..from.transect !=0 & !is.na(total_observations$Dist..from.transect), ]
total_observations <- total_observations$Dist..from.transect
total_observations <- cbind(c(1:length(total_observations)), rep(1, length(total_observations)), 
                            rep(1, length(total_observations)), total_observations)
total_observations <- data.frame(total_observations)
colnames(total_observations) <- c("object", "observer", "detected", "distance")


p_total <- ddf(method = 'ds',dsmodel =~ cds(key = "hr", formula=~1), 
               data = total_observations, meta.data = list(left = 50))
summary(p_total)
ddf.gof(p_total, main="Total observations goodness of fit")
plot(p_total)


#double observer for all species

distanceData <- function(species, overlap, lisa_missed, vic_missed) {
  
  #sets up data frame for distance sampling analysis
  #species code:  1 = baitfish, 2 = bottlenose dolphin
  
  dataframe <- cbind(1:nrow(overlap), rep(1, nrow(overlap)), rep(1, nrow(overlap)), overlap$Dist..from.transect, overlap$Species, overlap$Trial, overlap$Number)
  dataframe <- rbind(dataframe, dataframe)
  dataframe[(nrow(overlap) + 1):nrow(dataframe), 2] <- 2
  
  dataframe <- rbind(dataframe, cbind((nrow(dataframe) + 1):(nrow(dataframe)  + nrow(lisa_missed)), rep(2, nrow(lisa_missed)), rep(1, nrow(lisa_missed)), lisa_missed$Dist..from.transect, lisa_missed$Species, lisa_missed$Trial, lisa_missed$Number))
  dataframe <- rbind(dataframe, cbind((nrow(dataframe) + 1):(nrow(dataframe)  + nrow(lisa_missed)), rep(1, nrow(lisa_missed)), rep(0, nrow(lisa_missed)), lisa_missed$Dist..from.transect, lisa_missed$Species, lisa_missed$Trial, lisa_missed$Number))
  dataframe <- rbind(dataframe, cbind((nrow(dataframe) + 1):(nrow(dataframe)  + nrow(vic_missed)), rep(1, nrow(vic_missed)), rep(1, nrow(vic_missed)), vic_missed$Dist..from.transect, vic_missed$Species, vic_missed$Trial, vic_missed$Number))
  dataframe <- rbind(dataframe, cbind((nrow(dataframe) + 1):(nrow(dataframe)  + nrow(vic_missed)), rep(2, nrow(vic_missed)), rep(0, nrow(vic_missed)), vic_missed$Dist..from.transect, vic_missed$Species, vic_missed$Trial, vic_missed$Number))
  
  dataframe[, 1] <- c(rep(1:349, 2), rep(350:(350+nrow(lisa_missed) - 1), 2), rep((349 + 1 +nrow(lisa_missed)):(349 + nrow(lisa_missed) + nrow(vic_missed)), 2))
  
  dataframe <- data.frame(dataframe)
  colnames(dataframe) <- c("object", "observer", "detected", "distance", "species", "Trial", "size")
  
  dataframe <- dataframe[dataframe$species == species, ]
  
  if (species != 2) {
    dataframe <- dataframe[, 1:6]
  }
  
  dataframe <- dataframe[dataframe$distance !=0 & !is.na(dataframe$distance), ]

  return(dataframe)
  
}

calcAbundance <- function(area, dataframe, model) {
  
  #calculates density and abundance using distance sampling data
  #area = survey area (km2). Wollongong - Newcastle = 265km
  #dataframe = total_observations data frame from distanceData function
  #model = ddf model
  
  obs.table <- cbind(rep(1, nrow(dataframe)), dataframe$Trial, dataframe)
  colnames(obs.table)[1:2] <- c("Region.Label", "Sample.Label")
  
  region.table <- data.frame(matrix(c(1, area), ncol = 2, byrow = T))
  colnames(region.table) <- c("Region.Label", "Area")
  
  sample.table <- data.frame(cbind(rep(1, length(unique(dataframe$Trial))), unique(dataframe$Trial), rep(area, length(unique(dataframe$Trial)))))
  colnames(sample.table) <- c("Region.Label", "Sample.Label", "Effort")
  
  dht(model, region.table = region.table, sample.table = sample.table, obs.table = obs.table, 
      options = list(convert.units = 0.001))
  
}


total_observations <- distanceData(2, overlap, lisa_missed, vic_missed)

p_total <- ddf(method="io", mrmodel =~ glm(~distance), dsmodel =~ cds(key = "hr", formula=~1),
               data = total_observations, meta.data = list(left = 50, width = 500, point = FALSE))


summary(p_total)
ddf.gof(p_total, main="Total observations goodness of fit")
plot(p_total)


calcAbundance(265, total_observations, p_total)





