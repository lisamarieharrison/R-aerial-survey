#analysis of aerial survey data using only Lisa's observations
#author: Lisa-Marie Harrison
#date: 22/10/2015\

setwd("C:/Users/Lisa/Documents/phd/aerial survey/R/data")
lisa_obs <- read.csv("lisa_full_observations.csv", header = T)
library(mrds)
library(chron)
library(knitr)

table(lisa_obs$Species)

lisa_obs$Trial <- as.numeric(lisa_obs$Date)

#total sightings
table(lisa_obs$Flight.Direction, lisa_obs$Species)

#combine all shark species into a single category, except hammerheads
lisa_obs$Species[lisa_obs$Species %in% c("BS", "W", "Wh", "S")] <- "S"


#---------------------------- DETECTION FUNCTIONS -----------------------------#



createData <- function(species, lisa_obs, truncate=NULL, direction) {
  
  lisa_obs <- lisa_obs[lisa_obs$Flight.Direction == direction & lisa_obs$Type == "S" 
                       & lisa_obs$Dist..from.transect <= truncate &
                       lisa_obs$Dist..from.transect !=0 & !is.na(lisa_obs$Dist..from.transect)
                       & lisa_obs$Species == species, ]

  total_observations <- cbind(1:nrow(lisa_obs), rep(1, nrow(lisa_obs)), rep(1, nrow(lisa_obs)), lisa_obs$Dist..from.transect, lisa_obs$Species, lisa_obs$Trial, lisa_obs$Number)
  total_observations <- data.frame(total_observations)
  colnames(total_observations) <- c("object", "observer", "detected", "distance", "species", "Trial", "size")

  if (species != "BOT") {
    total_observations <- total_observations[, 1:6]
  } else {
    total_observations <- total_observations[!is.na(total_observations$size), ]
  }
  
  return(total_observations)
  
}

total_observations <- createData(species = "S", lisa_obs, truncate = 1000, direction = "S")

p_total <- ddf(method = 'ds',dsmodel =~ cds(key = "gamma", formula=~1), 
               data = total_observations, meta.data = list(left = 50, width = 1000))
summary(p_total)
ddf.gof(p_total, main="Total observations goodness of fit")
plot(p_total, main = "Baitfish - NORTH")



area <- 265*0.3
obs.table <- cbind(rep(1, nrow(total_observations)), total_observations$Trial, total_observations)
colnames(obs.table)[1:2] <- c("Region.Label", "Sample.Label")

region.table <- data.frame(matrix(c(1, area), ncol = 2, byrow = T))
colnames(region.table) <- c("Region.Label", "Area")

sample.table <- data.frame(cbind(rep(1, length(unique(total_observations$Trial))), unique(total_observations$Trial), rep(265, length(unique(total_observations$Trial)))))
colnames(sample.table) <- c("Region.Label", "Sample.Label", "Effort")

d <- dht(p_total, region.table = region.table, sample.table = sample.table, obs.table = obs.table, 
         options = list(convert.units = 0.001))
d





#-------------------------- Environmental variables ---------------------------#

#remove single "LT" or "RT" or times will be incorrect
test <- dat.south$Type[dat.south$Type %in% c("LT", "RT")]
for (i in 2:length(test)) {
  if ((test[i] == "LT" & test[i - 1] == "LT") | (test[i] == "RT" & test[i - 1] == "RT")) {
    dat.south <- dat.south[-i, ]
  }
}  

#initialize environmental variables
envt.var.south <- matrix(0, ncol = 27)
colnames(envt.var.south) <- c(
  "Beaufort.Sea.State.south.1",
  "Beaufort.Sea.State.south.2",
  "Beaufort.Sea.State.south.3",
  "Beaufort.Sea.State.south.4", 
  "Cloud.cover.south.0",
  "Cloud.cover.south.1",
  "Cloud.cover.south.2",
  "Cloud.cover.south.3",
  "Cloud.cover.south.4",
  "Cloud.cover.south.5",
  "Cloud.cover.south.6",
  "Cloud.cover.south.7",
  "Cloud.cover.south.8",
  "Water.clarity.south.1",
  "Water.clarity.south.2",
  "Water.clarity.south.3",
  "Glare.south.0",
  "Glare.south.10",
  "Glare.south.30",
  "Glare.south.15",
  "Glare.south.25",
  "Glare.south.5",
  "Glare.south.20",
  "Glare.south.35",
  "Glare.south.50",
  "Glare.south.60",
  "Glare.south.100"
)


dat.south <- lisa_obs[lisa_obs$Flight.Direction == "S", ]
dat.south <- dat.south[!dat.south$Time == "", ]

#sighting effort for southbound data
for (i in unique(dat.south$Date)) {
  for (j in 1:(nrow(dat.south[dat.south$Date == i, ]) - 1)) {
    
    if (dat.south[dat.south$Date == i, ]$Type[j] != "LT") {
      
      
      time.start <- chron(times. = as.character(dat.south$Time[dat.south$Date == i][j]), format = "h:m:s")
      time.stop  <- chron(times. = as.character(dat.south$Time[dat.south$Date == i][j + 1]), format = "h:m:s")
      mins <- hours(time.stop - time.start)*60 + minutes(time.stop - time.start)
      
      if (mins > 30) print(c(i, mins))
      
      #add minute differences to each environmental level
      for (k in c(22, 25, 26, 28)) {
        w <- which(colnames(envt.var.south) == as.name(paste(names(dat.south)[k], ".south.", dat.south[dat.south$Date == i, ][j, k], sep = "")))
        envt.var.south[1, w] <- envt.var.south[1, w] + mins
      }
    }
  }
}

sightingsAtLevel <- function(environmenal_hrs, dat, env_variable, species, flight_direction) {
  
  #calculates the number of sightings per hour at each environmental condition level
  #environmental_hrs = vector of hours at each level. e.g.: envt.var.south
  #dat = full data set of observations. e.g.: lisa_obs
  #env_variable = character containing the name of the environmental variable
  #species = character containing the species code. e.g.: "B" = baitfish
  #flight_direction = character containing the flight direction. "S" = south, "N" = north
  
  hrs <- environmenal_hrs[grep(env_variable, colnames(environmenal_hrs))]/60
  t <- table(factor(dat[dat$Species == species & dat$Flight.Direction == flight_direction, which(colnames(dat) == env_variable)],
                    levels = unique(dat[dat$Flight.Direction == flight_direction, 
                                        which(colnames(dat) == env_variable)])))
  t <- t[order(names(t))]                                                                                                                                 
  sightings_hr <- matrix(t/hrs, ncol = length(t))
  colnames(sightings_hr) <- names(t)
  
  #species hash table
  species_full <- cbind(c("B", "BOT", "S"), c("Baitfish", "Dolphins", "Sharks"))
  Species <- species_full[which(species_full[, 1] == species), 2]
  
  kable(sightings_hr, format = "pandoc", caption = paste("Number of sightings per hour effort at each", env_variable, "for", Species))
  
}

#sea state table
sightingsAtLevel(environmenal_hrs = envt.var.south, dat = lisa_obs, env_variable = "Beaufort.Sea.State", species = "B", flight_direction = "S")

#turbidity table
sightingsAtLevel(environmenal_hrs = envt.var.south, dat = lisa_obs, env_variable = "Water.clarity", species = "B", flight_direction = "S")


createTab <- function(species, dat, flight_direction, env_variable, env_dat=NULL) {
  
  #creates table for ggplot sightings per survey at environmental level
  #species = character with species code. e.g.: "B"
  #dat = full data frame of observations
  #flight_direction = character containing flight direction. "N" = north, "S" = south
  #env_variable = character containing environmental variable. e.g.: "Glare"
  #env_dat = vector of hours at each environmental level. Not needed if using Wind.Speed
  
  if (env_variable == "Wind.speed" | env_variable == "Wind.direction") {
    
    #wind separate because it is measured only once per survey
    
    wind_tab <- table(dat$Date[dat$Flight.Direction == flight_direction 
                               & dat$Species == species], 
                      factor(dat[dat$Flight.Direction == flight_direction 
                          & dat$Species == species, which(colnames(dat) == env_variable)], 
                          levels = unique(dat[dat$Flight.Direction == flight_direction,
                                              which(colnames(dat) == env_variable)])))
    wind_tab <- wind_tab[, order(as.numeric(colnames(wind_tab)))]
    wind_tab[wind_tab == 0] <- NA
    x <- as.numeric(colnames(wind_tab))
    y <- colMeans(wind_tab, na.rm = TRUE)
    y[is.nan(y)] <- 0
    
  } else {
    
    species_tab <- table(factor(dat[dat$Flight.Direction == flight_direction 
                                 & dat$Species == species, which(colnames(dat) == env_variable)], 
                             levels = unique(dat[dat$Flight.Direction == flight_direction,
                                                 which(colnames(dat) == env_variable)])))
    species_tab <- species_tab[order(names(species_tab))]
    cc_hrs <- env_dat[grep(env_variable, colnames(env_dat))]/60
    
    x <- as.numeric(names(species_tab))
    y <- c(species_tab/cc_hrs)
    
  }
  

  #create hash table for full species names
  species_full <- cbind(c("B", "BOT", "S"), c("Baitfish", "Dolphins", "Sharks"))
  
  Species <- rep(species_full[which(species_full[, 1] == species), 2], length(x))
  plot_dat <- data.frame(x, y, Species)
  
  return(plot_dat)
}


#wind speed
baitfish <- createTab(species = "B", dat = lisa_obs, flight_direction = "S", env_variable = "Wind.speed")
dolphins <- createTab(species = "BOT", dat = lisa_obs, flight_direction = "S", env_variable = "Wind.speed")
sharks   <- createTab(species = "S", dat = lisa_obs, flight_direction = "S", env_variable = "Wind.speed")

a <- rbind(baitfish, dolphins, sharks)
ggplot(a, aes(x = x,y = y, fill = Species, col = Species)) + 
  geom_point() + 
  geom_smooth(method = loess) + 
  xlab("Wind Speed (km/h)") + 
  ylab("Mean sightings per survey") + 
  scale_fill_manual(values = c("grey16", "red", "blue")) +
  scale_color_manual(values = c("grey16", "red", "blue")) + 
  theme(axis.text = element_text(colour = "black"), text = element_text(size = 30), legend.title = element_blank())


#Cloud cover
baitfish <- createTab(species = "B", dat = lisa_obs, flight_direction = "S", env_variable = "Cloud.cover", env_dat = envt.var.south)
dolphins <- createTab(species = "BOT", dat = lisa_obs, flight_direction = "S", env_variable = "Cloud.cover", env_dat = envt.var.south)
sharks   <- createTab(species = "S", dat = lisa_obs, flight_direction = "S", env_variable = "Cloud.cover", env_dat = envt.var.south)

a <- rbind(baitfish, dolphins, sharks)
ggplot(a, aes(x = x,y = y, fill = Species, col = Species)) + 
  geom_point() + 
  geom_smooth(method = loess) + 
  xlab("Cloud cover (0 = none, 8 = full)") + 
  ylab("Mean sightings per hour") + 
  scale_fill_manual(values = c("grey16", "red", "blue")) +
  scale_color_manual(values = c("grey16", "red", "blue")) + 
  theme(axis.text = element_text(colour = "black"), text = element_text(size = 30), legend.title = element_blank())






