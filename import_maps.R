#maps coastal NSW in R for use with aerial survey data 

#google earth maps using dismo library 
library(dismo)
g <- gmap(extent(150.8831, 151.7500, -34.4331, -32.9167), type = "satellite", lonlat = TRUE) 
plot(g)

#marmap library
library(marmap)
nsw <- getNOAA.bathy(lon1 = 150.8831,lon2 = 151.7500,lat1 = -34.4331,lat2 = -32.9167, resolution = 0.05)
plot(nsw, image = TRUE)

