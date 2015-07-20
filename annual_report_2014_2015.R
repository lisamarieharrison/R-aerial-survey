#annual 2014-2015 aerial survey report analysis
#author: Lisa-Marie Harrison
#date: 20/07/2015

setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/R/data")
dat <- read.csv("annual_interobserver_2014_2015.csv", header = T)
full_dat <- read.csv("annual_data_2014_2015.csv", header = T)
library(reshape)


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

#combine normal and interobserver sightings
dat <- rbind(full_dat, overlap, vic_missed)
dat$Number[is.na(dat$Number)] <- 1

#total sightings seen by each observer
table(dat$Species)

#table of sightings per season
table(dat$Flight.Direction[dat$Season == "Autumn"], dat$Species[dat$Season == "Autumn"])
table(dat$Flight.Direction[dat$Season == "Summer"], dat$Species[dat$Season == "Summer"])
table(dat$Flight.Direction[dat$Season == "Spring"], dat$Species[dat$Season == "Spring"])

#table of individuals per season
re <- dat[c(2, 7, 16, 17)]
names(re) <- c("Season", "Direction", "Species", "Number")
cast(re, Species ~ Season + Direction, sum)

#remove species of no interest
dat <- dat[!(dat$Species %in% c("J", "HB", "R", "T", "COM", "SF", "F", "P", "")), ]
dat$Species <- factor(dat$Species)
dat$Number[is.na(dat$Number)] <- 1

#sightings per day
days_per_season <- c(8, 8, 12)
sightings_per_day <- round(table(dat$Season, dat$Species)/days_per_season, 2)

par(mar = c(10, 5, 2, 2))
barplot(sightings_per_day, xaxt = "n", beside = T, col = c("gray20", "grey60", "grey85"), ylab = "Mean sightings per day")
legend("topright", c("autumn (n = 8)", "spring (n = 8)", "summer (n = 12)"), fill = c("gray20", "grey60", "grey85"), bty = "n")
axis(1, at = c(2, 6, 10, 14, 18, 22, 26), las = 2, labels= c("Baitfish", "Bottlenose dolphin", "Hammerhead shark", "Unid. shark", "White shark", "Bull shark", "Whaler shark"))


#individuals per sighting
re <- dat[c(2, 16, 17)]
names(re) <- c("Season", "Species", "Number")
num_sightings <- cast(re, Species ~ Season)[, 2:4]
num_inds <- cast(re, Species ~ Season, sum)[, 2:4]
inds_per_sighting <- num_inds/num_sightings
inds_per_sighting_tab <- as.table(t(inds_per_sighting))
colnames(inds_per_sighting_tab) <- cast(re, Species ~ Season)[, 1]

par(mar = c(10, 5, 2, 2))
barplot(inds_per_sighting_tab, xaxt = "n", beside = T, col = c("gray20", "grey60", "grey85"), ylab = "Mean sightings per day", legend = T, args.legend = c(bty= "n"))
axis(1, at = c(2, 6, 10, 14, 18, 22, 26), las = 2, labels= c("Baitfish", "Bottlenose dolphin", "Hammerhead shark", "Unid. shark", "White shark", "Bull shark", "Whaler shark"))


