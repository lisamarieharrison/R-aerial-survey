#test distance sampling assumptions using both directions
#author: Lisa-Marie Harrison
#date: 13/12/2016

if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/aerial survey/R")
} else {
  setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/R")
}


lisa_obs <- read.csv("data/lisa_full_observations.csv", header = T)
library(mrds)
library(chron)
library(knitr)
library(ggplot2)
library(reshape2)
library(gridExtra)

#add a unique trial number to each day
date_levels <- as.numeric(as.factor(c(as.character(lisa_obs$Date))))
lisa_obs$Trial <- date_levels[1:nrow(lisa_obs)]

lisa_obs$Year <- substr(lisa_obs$Date, nchar(as.character(lisa_obs$Date)) - 3, nchar(as.character(lisa_obs$Date)))

#source other required functions
file_list <- c("mrds_modified_functions.R", 
               "calcAbundanceSingle.R",
               "sightingsAtLevel.R",
               "createDistanceData.R",
               "corObsByPercent.R",
               "calcEnvtEffort.R",
               "correctLTRT.R")

for (f in file_list) {
  source(paste("code/R-aerial-survey/functions/", f, sep =""))
}

#merge lisa and vic's observations to get all sightings
table(lisa_obs$Flight.Direction[!(lisa_obs$Species == "")], lisa_obs$Species[!(lisa_obs$Species == "")])

poolAllSharks <- function(x) {
  
  x$Species[x$Species %in% c("W", "S", "Wh", "BS")] <- "S"
  
  return(x)
  
}

#combine all shark species into a single category, except hammerheads
lisa_obs <- poolAllSharks(lisa_obs)


#---------------------------- DETECTION FUNCTIONS -----------------------------#

strip_width <- 300


#plots of detection probability for each species (300m strip width)
for (s in c("B", "BOT", "S")) {
  par(mfrow = c(1, 2))
  
  for (d in c("N", "S")) {
    
    total_observations <- createDistanceData(species = s, lisa_obs, truncate = strip_width, direction = d)
    
    p_total <- ddf(method = 'ds',dsmodel =~ cds(key = "gamma", formula=~1), 
                   data = total_observations, meta.data = list(left = 50, width = 300))
    summary(p_total)
    
    #species and flight direction hash tables
    species_tab <- matrix(c("B", "BOT", "S", "Baitfish", "Bottlenose dolphins", "Sharks"), ncol = 2)
    dir_tab <- matrix(c("N", "S", "North", "South"), ncol = 2)
    
    plot(p_total, main = , cex.axis = 2, cex.lab = 2)
    title(paste(species_tab[which(species_tab[, 1] == s), 2], ": ", dir_tab[which(dir_tab[, 1] == d), 2], sep = ""), 
          cex.main = 3)
    
  }
}
