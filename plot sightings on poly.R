setwd("C:/Users/43439535/Documents/Lisa/phd/aerial survey/data/speed trials")

lisa  <- read.csv("speed_trials_datasheet_lisa.csv", header = T, skip = 1)
rob   <- read.csv("speed_trials_datasheet_rob.csv", header = T, skip = 1)
vic   <- read.csv("speed_trials_datasheet_vic.csv", header = T, skip = 1)
sally <- read.csv("speed_trials_datasheet_sally.csv", header = T, skip = 1)

library(zoom)
library(mapdata)
library(ggplot2)
library(plotly)
library(maps)
library(maptools)

nsw  <- map("world2Hires", regions=c("Australia"), ylim = c(-32.9167, -32.3333), xlim = c(151.7500, 152.5333), plot = FALSE)
nsw_poly <- map2SpatialPolygons(nsw, IDs = c("1", "2"), proj4string=CRS("+proj=longlat +datum=WGS84"))
nsw_ggplot   <- fortify(nsw_poly, region="id") #df for ggplot


trial <- 1

lisa_temp <- subset(lisa, lisa$Trial_number == trial & lisa$Type == "S")
vic_temp <- subset(vic, vic$Trial_number == trial & vic$Type == "S")
rob_temp <- subset(rob, rob$Trial_number == trial & rob$Type == "S")
sally_temp <- subset(sally, sally$Trial_number == trial & sally$Type == "S")

lisa_points <- cbind(lisa_temp$Long, lisa_temp$Lat)
lisa_labels <- paste(lisa_temp$Species, lisa_temp$Number)

vic_points  <- cbind(vic_temp$Long, vic_temp$Lat)
vic_labels  <- paste(vic_temp$Species, vic_temp$Number)

rob_points  <- cbind(rob_temp$Long, rob_temp$Lat)
rob_labels  <- paste(rob_temp$Species, rob_temp$Number)

sally_points  <- cbind(sally_temp$Long, sally_temp$Lat)
sally_labels  <- paste(sally_temp$Species, sally_temp$Number)


p <- ggplot() + 
  geom_polygon(data=nsw_ggplot, aes(x=long, y=lat, group=id), color="black", fill = "grey") +
  xlim(151.7500, 152.5333) + 
  ylim(-32.9167, -32.3333) +
  geom_text(aes(lisa_points[, 1], lisa_points[, 2]), label = lisa_labels, colour = "red") + 
  geom_text(aes(vic_points[, 1], vic_points[, 2]), label = vic_labels, colour = "blue") +
  geom_text(aes(rob_points[, 1], rob_points[, 2]), label = rob_labels, colour = "green") +
  geom_text(aes(sally_points[, 1], sally_points[, 2]), label = sally_labels, colour = "orange")
ggplotly(p)



