# plot shark attack locations around australia
# author: Lisa-Marie Harrison

setwd("C:/Users/43439535/Documents/Lisa/phd/aerial survey/data/speed trials")

dat  <- read.csv("speed_trials_datasheet_lisa.csv", header = T, skip = 1)

#install plotly using devtools
install.packages("viridis") # dependency
install.packages("devtools")
devtools::install_github("ropensci/plotly")

library(mapdata)
library(maps)
library(maptools)
library(ggplot2)
library(plotly)

aust  <- map("world2Hires", regions="Australia", fill = TRUE, col = "grey", xlim = c(95, 160), ylim = c(-45, -10.05140))
aust_poly <- map2SpatialPolygons(aust, IDs = aust$names, proj4string=CRS("+proj=longlat +datum=WGS84"))
aust_ggplot   <- fortify(aust_poly, region="id") #df for ggplot


p <- ggplot() + 
  geom_polygon(data=aust_ggplot, aes(x=long, y=lat, group=group), color="black", fill = "grey") 
  #geom_point(aes(dat$longitude, dat$latitude, shape = dat$species, colour = dat$species)) 
ggplotly(p)
