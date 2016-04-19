#analysis of aerial survey data using only Lisa's observations
#author: Lisa-Marie Harrison
#date: 22/10/2015

if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/aerial survey/R")
} else {
  setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/R")
}


lisa_obs <- read.csv("data/lisa_full_observations.csv", header = T)
vic_obs  <- read.csv("data/vic_extra_observations.csv", header = T)
library(mrds)
library(chron)
library(knitr)
library(ggplot2)
library(reshape2)
library(gridExtra)

#add a unique trial number to each day
date_levels <- as.numeric(as.factor(c(as.character(lisa_obs$Date), as.character(vic_obs$Date))))
lisa_obs$Trial <- date_levels[1:nrow(lisa_obs)]
vic_obs$Trial <- date_levels[(1 + nrow(lisa_obs)):(nrow(vic_obs) + (nrow(lisa_obs)))]

#source other required functions
file_list <- c("mrds_modified_functions.R", 
               "calcAbundanceSingle.R",
               "sightingsAtLevel.R",
               "createDistanceData.R",
               "corObsByPercent.R",
               "calcEnvtEffort.R",
               "correctLTRT.R")

for (f in file_list) {
  source(paste("R-aerial-survey/functions/", f, sep =""))
}

#merge lisa and vic's observations to get all sightings
all_obs <- rbind(lisa_obs, vic_obs)
table(all_obs$Flight.Direction, all_obs$Species)

season <- rbind(lisa_obs, vic_obs)
annual <- season

poolAllSharks <- function(x) {
  
  x$Species[x$Species %in% c("W", "S", "Wh", "BS")] <- "S"
  
  return(x)
  
}

#combine all shark species into a single category, except hammerheads
lisa_obs <- poolAllSharks(lisa_obs)
vic_obs  <- poolAllSharks(vic_obs)

dat.south <- rbind(lisa_obs[lisa_obs$Flight.Direction == "S", ], vic_obs[vic_obs$Flight.Direction == "S", ])

#correct vics number of sightings using percentage of obs seen by lisa
source("~/Lisa/phd/aerial survey/R/R-aerial-survey/percentOfLisa.R")
percent_of_lisa <- percentOfLisa()
vic_obs <- corObsByPercent(dat = vic_obs, percent_mat = percent_of_lisa) 

#combine lisa and vics observations to get all corrected observations
cor_obs <- rbind(lisa_obs, vic_obs)

#---------------------------- DETECTION FUNCTIONS -----------------------------#

strip_width <- 1000


#plots of detection probability for each species (300m strip width)
for (s in c("B", "BOT", "S")) {
  par(mfrow = c(1, 2))
  
  for (d in c("N", "S")) {
    
    total_observations <- createDistanceData(species = s, cor_obs, truncate = strip_width, direction = d)
    
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

#abundance calculator

for (s in c("B", "BOT", "S")) {
  
  total_observations <- createDistanceData(species = s, cor_obs, truncate = strip_width, direction = "S")
  
  p_total <- ddf(method = 'ds',dsmodel =~ cds(key = "gamma", formula=~1), 
                 data = total_observations, meta.data = list(left = 50, width = strip_width))
  print(summary(p_total))
  #ddf.gof(p_total, main="Total observations goodness of fit")
  
  d <- calcAbundanceSingle(265*strip_width/1000, total_observations, p_total)
  print(d)
  
}

#-------------------------- Environmental variables ---------------------------#

#remove single "LT" or "RT" or times will be incorrect
dat.south <- correctLTRT(dat.south)


envt.var.south <- calcEnvtEffort(dat.south)

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
surveys_per_season <- colSums(table(season$Trial, substr(season$Season, 1, 6)) > 0)*2
round(table(substr(season$Season, 1, 6), season$Species)/surveys_per_season, 1)

#combine all shark species into one category for plotting
season$Species[season$Species %in% c("BS", "W", "Wh", "S")] <- "S"

#boxplots of main species by season
for (s in c("B", "BOT", "S")) {
  
  season_tab <- table(substr(season$Season, 1, 6)[season$Species == s], season$Trial[season$Species == s])
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


#----------------------------ANNUAL DIFFERENCES--------------------------------#

surveys_per_season <- colSums(table(season$Trial, season$Season) > 0)*2
round(table(season$Season, season$Species)/surveys_per_season, 1)

for (s in c("B", "BOT", "S")) { 
  
  season_tab <- table(cor_obs$Season[cor_obs$Species == s], cor_obs$Trial[cor_obs$Species == s])
  season_tab[season_tab == 0] <- NA
  
  #species hash table
  species_hash <- matrix(c("B", "BOT", "S", "Baitfish", "Bottlenose dolphins", "Sharks"), ncol = 2)
  
  print(ggplot(data = melt(t(season_tab)), aes(x=Var2, y=value)) + geom_boxplot(aes(fill=Var2), alpha = 0.5) + 
    xlab("") + 
    ylab("Sightings per survey") + 
    scale_fill_manual(name = "Season", values = c("red", "blue", "yellow")) + 
    ggtitle(paste(species_hash[which(species_hash[, 1] == s), 2])) +
    theme(plot.title = element_text(lineheight = 1, face="bold", size = 30), axis.text = element_text(size = 25),
          legend.position="none", axis.title = element_text(size = 25)))

}



