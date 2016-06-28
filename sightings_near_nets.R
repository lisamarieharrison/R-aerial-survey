#sightings near nets
#author: Lisa-Marie Harrison
#date: 27/06/2016

if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/aerial survey/data")
  source_location <- "~/Lisa/phd/Mixed models/R code/R-functions-southern-ocean/"
} else {
  setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/data")
  source_location <- "~/phd/southern ocean/Mixed models/R code/R-functions-southern-ocean/"
}

dat  <- read.csv("aerial_survey_summary_r.csv", header = T)
nets <- read.csv("shark_net_locations.csv", header = T)

source_list <- c("deg2rad.R",
       "gcdHF.R")

for (f in source_list) {
  source(paste0(source_location, f))
}

dat <- dat[dat$Type == "S", ]
dat$Lat <- as.numeric(as.character(dat$Lat))
dat <- dat[!is.na(dat$Lat), ]
dat <- dat[dat$Species == "BOT", ]

animals <- NULL
area    <- NULL


inclusionPercent <- function(rad) {
  
  #how to account for overlap of radiuses when net radiuses are overlapping?
  
  lats_to_check <- seq(min(nets$latitude), max(nets$latitude), by = 0.001) #check 100m sequence
  point_in_nets <- rep(FALSE, length(lats_to_check))
  for (i in 1:length(lats_to_check)) {
    
    for (j in 1:nrow(nets)) {
      
      dist_to_net <- gcdHF(deg2rad(lats_to_check[i]), deg2rad(nets$longitude[j]), deg2rad(nets$latitude[j]), deg2rad(nets$longitude[j]))*1000 #m
      
      if (dist_to_net <= rad) {
        
        point_in_nets[i] <- TRUE
        
        break
        
      }
      
      
    }  
    
  }
  
  percentage_coverage <- sum(point_in_nets)/length(point_in_nets)*100
  
  return (percentage_coverage)
  
}

for (radius in seq(100, 4000, by = 100)) {
  
  radius_around_nets <- radius #m
  
  at_beach <- matrix(FALSE, nrow = nrow(nets), ncol = nrow(dat))
  for (i in 1:ncol(at_beach)) {
    
    for (j in 1:nrow(at_beach)) {
      
      distance <- gcdHF(deg2rad(dat$Lat[i]), deg2rad(nets$longitude[j]), deg2rad(nets$latitude[j]), deg2rad(nets$longitude[j]))*1000 #m
      
      if (distance < radius_around_nets) {
        
        at_beach[j, i] <- TRUE
        
        break
        
      }
      
    }
    
  }
  
  #percent of sightings within specified radius of nets
  animals <- c(animals, table(colSums(at_beach))[2] / sum( table(colSums(at_beach)))*100)
  
  #survey length = 265km
  #covered area by nets 
  
  area <- c(area, inclusionPercent(radius_around_nets))
  #area <- c(area, (51*2*radius_around_nets/1000)/265 * 100)
  
}

plot(area, animals, xlim = c(1, 100), ylim = c(1, 100))
points(c(0, 100), c(0, 100), type = "l", col = "red")





#closest nets to each other to check how soon overlap will occur
net_distance <- matrix(NA, nrow = nrow(nets), ncol = nrow(nets))
for (j in 1:nrow(nets)) {
  
  for (i in 1:nrow(nets)) {
    
    net_distance[j, i] <- gcdHF(deg2rad(nets$latitude[i]), deg2rad(nets$longitude[i]), deg2rad(nets$latitude[j]), deg2rad(nets$longitude[j]))*1000 #m
    
  }
  
}
diag(net_distance) <- NA
apply(net_distance, 2, min, na.rm = TRUE) #some nets have other nets only 180m away
