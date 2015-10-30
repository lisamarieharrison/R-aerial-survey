#analysis of aerial survey data using only Lisa's observations
#author: Lisa-Marie Harrison
#date: 22/10/2015

setwd("C:/Users/Lisa/Documents/phd/aerial survey/R/data")
lisa_obs <- read.csv("lisa_full_observations.csv", header = T)
vic_obs  <- read.csv("vic_extra_observations.csv", header = T)
library(mrds)
library(chron)
library(knitr)
library(ggplot2)
library(reshape2)
library(gridExtra)

#source mrds code modified to use only single sided strip width
source("C:/Users/Lisa/Documents/phd/aerial survey/R/code/R-aerial-survey/mrds_modified_functions.R")

#merge lisa and vic's observations to get all sightings
all_obs <- rbind(lisa_obs, vic_obs)
table(all_obs$Flight.Direction, all_obs$Species)

season <- rbind(lisa_obs, vic_obs)

#combine all shark species into a single category, except hammerheads
lisa_obs$Species[lisa_obs$Species %in% c("BS", "W", "Wh", "S")] <- "S"
vic_obs$Species[vic_obs$Species %in% c("BS", "W", "Wh", "S")] <- "S"

dat.south <- rbind(lisa_obs[lisa_obs$Flight.Direction == "S", ], vic_obs[vic_obs$Flight.Direction == "S", ])

for (t in unique(vic_obs$Trial)) {
  for (s in c("B", "BOT", "S")) {
    cor_factor <- percent_of_lisa[which(names(percent_of_lisa) == s)]
    n_rows <- nrow(vic_obs[vic_obs$Trial == t & vic_obs$Species == s, ])
    if (n_rows != 0) {
      n_rep <- abs(n_rows - round(n_rows / cor_factor))
      rep_row <- vic_obs[vic_obs$Trial == t & vic_obs$Species == s, ][1, ]
      rep_row$Dist..from.transect <- mean(na.omit(vic_obs[vic_obs$Trial == t & vic_obs$Species == s, ]$Dist..from.transect))
      for (i in 1:n_rep) {
        vic_obs <- rbind(vic_obs, rep_row)
      }
    }
  }
}

cor_obs <- rbind(lisa_obs, vic_obs)

#---------------------------- DETECTION FUNCTIONS -----------------------------#



createData <- function(species, lisa_obs, direction, truncate=NULL) {
  
  #creates data frame in the form mrds requires for distance sampling
  #species = character with the species code. e.g.: "BOT" = bottlenose dolphin
  #lisa_obs = data.frame of all observations
  #direction = character with flight direction. "N" = north, "S" = south
  #truncate = optional numeric specifying the distance (m) from transect to truncate sightings
  #return = data.frame with columns object, observer, detected, distance, species, Trial, size
  
  lisa_obs <- lisa_obs[lisa_obs$Flight.Direction == direction & 
                         lisa_obs$Type == "S" & 
                         lisa_obs$Dist..from.transect !=0 & 
                         !is.na(lisa_obs$Dist..from.transect) & 
                         lisa_obs$Species == species, ]

  if (!is.null(truncate)) {
    lisa_obs <- lisa_obs[lisa_obs$Dist..from.transect <= truncate, ]
  }
  
  total_observations <- cbind(1:nrow(lisa_obs), rep(1, nrow(lisa_obs)), rep(1, nrow(lisa_obs)), lisa_obs$Dist..from.transect, lisa_obs$Species, lisa_obs$Trial, lisa_obs$Number)
  total_observations <- data.frame(total_observations)
  colnames(total_observations) <- c("object", "observer", "detected", "distance", "species", "Trial", "size")

  total_observations <- apply(total_observations, 2, as.character)
  total_observations <- data.frame(apply(total_observations, 2, as.numeric))
  
  if (species == "BOT") {
    total_observations <- total_observations[!is.na(total_observations$size), ]
  } else {
    total_observations <- total_observations[, 1:6]
  }
  
  return(total_observations)
  
}

calcAbundanceSingle <- function(area, dataframe, model) {
  
  #calculates density and abundance using distance sampling data
  #area = survey area (km2). Wollongong - Newcastle = 265km
  #dataframe = total_observations data frame from distanceData function
  #model = ddf model
  
  obs.table <- cbind(rep(1, nrow(dataframe)), dataframe$Trial, dataframe)
  colnames(obs.table)[1:2] <- c("Region.Label", "Sample.Label")
  
  region.table <- data.frame(matrix(c(1, area), ncol = 2, byrow = T))
  colnames(region.table) <- c("Region.Label", "Area")
  
  sample.table <- data.frame(cbind(rep(1, 54), c(1:54), rep(265, 54)))
  colnames(sample.table) <- c("Region.Label", "Sample.Label", "Effort")
  
  d <- dht(model, region.table = region.table, sample.table = sample.table, obs.table = obs.table, 
           options = list(convert.units = 0.001))
  
  return(d)
  
}

strip_width <- 1000

total_observations <- createData(species = "S", cor_obs, truncate = strip_width, direction = "S")

p_total <- ddf(method = 'ds',dsmodel =~ cds(key = "gamma", formula=~1), 
               data = total_observations, meta.data = list(left = 50, width = strip_width))
summary(p_total)
#ddf.gof(p_total, main="Total observations goodness of fit")

d <- calcAbundanceSingle(265*strip_width/1000, total_observations, p_total)
d


#plots of detection probability for each species (300m strip width)
for (s in c("B", "BOT", "S")) {
  par(mfrow = c(1, 2))
  
  for (d in c("N", "S")) {
    
    total_observations <- createData(species = s, cor_obs, truncate = strip_width, direction = d)
    
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




#-------------------------- Environmental variables ---------------------------#

#remove single "LT" or "RT" or times will be incorrect
dat.south <- dat.south[!dat.south$Time == "", ]
dat.south[, c("Beaufort.Sea.State", "Water.clarity")] <- apply(dat.south[, c("Beaufort.Sea.State", "Water.clarity")], 2, as.integer)
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

#sighting effort for southbound data
for (i in unique(dat.south$Date)) {
  for (j in 1:(nrow(dat.south[dat.south$Date == i, ]) - 1)) {
    
    if (dat.south[dat.south$Date == i, ]$Type[j] != "LT") {
      
      
      time.start <- chron(times. = as.character(dat.south$Time[dat.south$Date == i][j]), format = "h:m:s")
      time.stop  <- chron(times. = as.character(dat.south$Time[dat.south$Date == i][j + 1]), format = "h:m:s")
      mins       <- hours(time.stop - time.start)*60 + minutes(time.stop - time.start)
      
      if (mins > 30) {
        print(c(i, mins)) #check for times that have been typed in incorrectly
      }
      
      #add minute differences to each environmental level
      for (k in c(24, 27, 28, 30)) {
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
sightingsAtLevel(environmenal_hrs = envt.var.south, dat = dat.south, env_variable = "Beaufort.Sea.State", species = "S", flight_direction = "S")

#turbidity table
sightingsAtLevel(environmenal_hrs = envt.var.south, dat = dat.south, env_variable = "Water.clarity", species = "S", flight_direction = "S")


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
  coord_cartesian(ylim = c(0, 45)) +
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
  coord_cartesian(ylim = c(0, 45)) +
  scale_fill_manual(values = c("grey16", "red", "blue")) +
  scale_color_manual(values = c("grey16", "red", "blue")) + 
  theme(axis.text = element_text(colour = "black"), text = element_text(size = 30), legend.title = element_blank())


#------------------------PERCEPTION BIAS FOR SHARKS----------------------------#

#Robbins et al found 17.1% of shark analogues were seen from a helicopter

est <- 2.4 #uncorrected abundance estimate for 1000m S transect
adjusted <- est/0.171


#-----------------------AVAILABILITY BIAS FOR DOLPHINS-------------------------#

#Gazo et al 2004 found that mean proportion of time that a bottlenose dolphin pod
#spent at the surface was 0.77. 23% of pods are underwater and unavailable for sampling
#173 groups seen with average group size = 16.1

adjusted <- (173/0.77)/(265*54*0.76*0.32)*16.1*265


#--------------------------SEASONAL DIFFERENCES--------------------------------#

#number of surveys per season
surveys_per_season <- colSums(table(season$Trial, season$Season) > 0)*2
round(table(season$Season, season$Species)/surveys_per_season, 1)

#boxplots of main species by season
for (s in c("B", "BOT", "S")) {
  
  season_tab <- table(season$Season[season$Species == s], season$Trial[season$Species == s])
  season_tab[season_tab == 0] <- NA
  
  #species hash table
  species_hash <- matrix(c("B", "BOT", "S", "Baitfish", "Bottlenose dolphins", "Sharks", "Sightings per survey", "", ""), ncol = 3)
  
  assign(paste("p", which(species_hash[, 1] == s), sep=""), 
    ggplot(data = melt(t(season_tab)), aes(x=Var2, y=value)) + geom_boxplot(aes(fill=Var2), alpha = 0.5) + 
    xlab("") + 
    ylab(species_hash[which(species_hash[, 1] == s), 3]) + 
    scale_y_continuous(limits = c(0, 80)) +
    scale_fill_manual(name = "Season", values = c("red", "blue", "yellow")) + 
    ggtitle(paste(species_hash[which(species_hash[, 1] == s), 2])) +
    theme(plot.title = element_text(lineheight = 1, face="bold", size = 30), axis.text = element_text(size = 25),
          legend.position="none", axis.title = element_text(size = 25)))
  
}

plot_list <- list(p1, p2, p3)
do.call(grid.arrange, c(plot_list, list(ncol = 3))) #plot 3 plots next to each other







