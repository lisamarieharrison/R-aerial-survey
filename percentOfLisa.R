percentOfLisa <- function() {
  
  #calculate the percentage of lisa's observations that Vic saw
  
  dat <- read.csv("~/Lisa/phd/aerial survey/R/data/interobserver_20150407.csv", header = T)
  
  #source required functions
  source("~/Lisa/phd/aerial survey/R/R-aerial-survey/functions/checkSameSighting.R")
  source("~/Lisa/phd/aerial survey/R/R-aerial-survey/functions/findOverlappingObservations.R")
  
  #dat <- dat[dat$Type == "S", ]
  dat <- dat[dat$Species != "", ]
  dat$Species[dat$Species == "B "] <- "B"
  dat$Species <- as.factor(dat$Species)
  
  lisa <- dat[dat$Observer == "Lisa", ]
  vic  <- dat[dat$Observer == "Vic", ]
  
  #separate observations into overlap (both seen) and missed observations
  obs <- findOverlappingObservations(lisa, vic)
  overlap     <- obs$overlap
  lisa_missed <- obs$obs1_missed
  vic_missed  <- obs$obs2_missed
  
  #total number of sightings by each observer
  lisa_missed <- poolAllSharks(lisa_missed)
  overlap     <- poolAllSharks(overlap)
  vic_missed  <- poolAllSharks(vic_missed)
  
  lisa_missed <- lisa_missed[lisa_missed$Species %in% levels(overlap$Species), ]

  vic_seen  <- table(lisa_missed$Species) + table(overlap$Species)
  lisa_seen <- table(vic_missed$Species) + table(overlap$Species)
  
  #vic percentage of Lisa
  percent_of_lisa <- vic_seen/lisa_seen
  
  return(percent_of_lisa)
  
}
