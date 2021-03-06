#dolphin sightings around shark nets
#author: Lisa-Marie Harrison
#date: 27/06/2016

if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/aerial survey/data")
  source_location <- "~/Lisa/phd/Mixed models/R code/R-functions-southern-ocean/"
} else {
  setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/data")
  source_location <- "~/phd/southern ocean/Mixed models/R code/R-functions-southern-ocean/"
}
library(sp)
library(spatstat)
library(plyr)

dat     <- read.csv("aerial_survey_summary_r.csv", header = T) #lisa's sighting data
nets    <- read.csv("shark_net_locations.csv", header = T) #coordinates of shark nets
no_nets <- read.csv("beaches_without_nets.csv", header = T) #coordinates of sandy beaches with no shark net

source_list <- c("deg2rad.R",
                 "gcdHF.R")

invisible(Map(source, paste0(source_location, source_list)))


dat <- dat[dat$Type == "S", ]
dat$Lat  <- as.numeric(as.character(dat$Lat))
dat$Long <- as.numeric(as.character(dat$Long))
dat <- dat[!is.na(dat$Lat), ]
dat$Lat[dat$Lat > 0] <- dat$Lat[dat$Lat > 0]*-1
dat$Transect <- as.numeric(dat$Date)

dat_bot  <- dat[dat$Species == "BOT", ]
dat_fish <- dat[dat$Species == "B", ]
dat_fish <- dat_fish[dat_fish$Length.Size != "S", ]

distAllCombinations <- function (lat1, lat2, long) {
  
  coords_to_compare <- deg2rad(data.frame("lat1" = rep(lat1, each = length(lat2)), "long1" = rep(mean(long), each = length(lat2)), "lat2" = lat2, "long2" = mean(long)))
  dist_all <- mapply(gcdHF, coords_to_compare$lat1, coords_to_compare$long1, coords_to_compare$lat2, coords_to_compare$long2)*1000 #m
  
  dist_all_mat <- matrix(dist_all, nrow = length(lat1), byrow = TRUE) 
  
  return (dist_all_mat)
  
}


withinRadiusOfBeach <- function(radius, net_coords, lats, longs=net_coords$longitude) {
  
  #function to check percentage coverage of survey area by shark net_coords
  #radius = vector of distances (m) of radius from shark net for a sighting to be considered to be at the net
  #lats = vector of latitudes to check if they are near the net_coords
  #longs = optional vector of longitudes. Defaults to using net_coords if not specified
  #return: numeric of percentage of points that fell within the specified distance to the net

  point_in_net_coords <- rep(FALSE, length(lats)) #boolean: is the point near a net
  
  #coords_to_compare <- deg2rad(data.frame("lat1" = rep(lats, each = nrow(net_coords)), "long1" = longs, "lat2" = net_coords$latitude, "long2" = longs))
  #dist_to_net <- mapply(gcdHF, coords_to_compare$lat1, coords_to_compare$long1, coords_to_compare$lat2, coords_to_compare$long2)*1000 #m
  
  dist_to_net_mat <- distAllCombinations(lats, net_coords$latitude, longs)
  
  #dist_to_net_mat <- matrix(dist_to_net, nrow = length(lats), byrow = TRUE) 
  
  point_at_net <-lapply(radius, function(x) rowSums(dist_to_net_mat <= x) > 0)

  perc_in_radius <- unlist(lapply(point_at_net, function(x)  sum(x)/length(x)*100))
  
  return (perc_in_radius)
  
}

#radiuses to check
radius <- seq(100, 4000, by = 100)

#percentage of animals at beach
animals <- withinRadiusOfBeach(radius, nets, dat_bot$Lat)

#covered area by nets at each radius
lats_to_check <- seq(min(nets$latitude), max(nets$latitude), by = 0.001) #check 100m sequence
area <- withinRadiusOfBeach(radius, nets, lats_to_check)

plot(area, animals, xlim = c(1, 100), ylim = c(1, 100))
points(c(0, 100), c(0, 100), type = "l", col = "red")



#closest nets to each other to check how soon overlap will occur
net_distance  <- distAllCombinations(nets$latitude, nets$latitude, nets$longitude)
diag(net_distance) <- NA
apply(net_distance, 2, min, na.rm = TRUE) #some nets have other nets only 180m away


#distance of each sighting to closest net
sighting_to_net  <- distAllCombinations(dat_bot$Lat, nets$latitude, nets$longitude)
#hist(apply(sighting_to_net, 2, min, na.rm = TRUE), xlab = "min distance to net (m)")


# ----------------------- SPATIAL ANALYSIS OF DOLPHINS AND NETS -------------------------#

marked_points <- data.frame("longitude" = dat_bot$Long, "latitude" = dat_bot$Lat, "mark" = "dolphin", "seg" = dat_bot$Transect)
marked_points <- rbind(marked_points, data.frame("longitude" = dat_fish$Long, "latitude" = dat_fish$Lat, "mark" = "fish", "seg" = dat_fish$Transect))
marked_points <- na.omit(marked_points)

#convert to utm so units are in m
xy = data.frame(marked_points$longitude, marked_points$latitude)
colnames(coordinates(xy)) <- c("lon", "lat")
proj4string(xy) <- CRS("+proj=longlat +dat_botum=WGS84")
marked_points[, 1:2] <- coordinates(spTransform(xy, CRS("+proj=utm +zone=53 ellps=WGS84")))


#Single line segment for each transect

marked_points$tp <- (marked_points$latitude - min(marked_points$latitude))/(max(marked_points$latitude) - min(marked_points$latitude))
owin <- owin(xrange = c(mean(marked_points$longitude), mean(marked_points$longitude)+500), yrange = c(min(marked_points$latitude), max(marked_points$latitude)))
transect_bounds <- ppp(x = rep(mean(marked_points$longitude), 47*2), y = rep(c(min(marked_points$latitude), max(marked_points$latitude)), 47), window = owin)
jointed_vertices <- matrix(FALSE, nrow = 94, ncol = 94)
pairs <- cbind(seq(1, 93), seq(2, 94))
jointed_vertices[rbind(pairs, cbind(pairs[, 2], pairs[, 1]))] <- TRUE


transect_line <- linnet(vertices = transect_bounds, m = jointed_vertices)
marked_points$longitude <- mean(transect_bounds$x)
names(marked_points) <- c("x", "y", "mark", "seg", "tp")
point_network <- lpp(X = marked_points, L = transect_line)


#dolphin fish school correlation
pair_correlation <- linearpcfcross(point_network, "fish", "dolphin")
plot(pair_correlation) 

K <- linearKcross(point_network, "fish", "dolphin")
plot(K)


#single day

r <- NA
est <- NA

for (day in unique(marked_points$seg)) {
  
  marked_points <- data.frame("longitude" = dat_bot$Long, "latitude" = dat_bot$Lat, "mark" = "dolphin", "seg" = dat_bot$Transect)
  marked_points <- rbind(marked_points, data.frame("longitude" = dat_fish$Long, "latitude" = dat_fish$Lat, "mark" = "fish", "seg" = dat_fish$Transect))
  marked_points <- na.omit(marked_points)

  #convert to utm so units are in m
  xy = data.frame(marked_points$longitude, marked_points$latitude)
  colnames(coordinates(xy)) <- c("lon", "lat")
  proj4string(xy) <- CRS("+proj=longlat +dat_botum=WGS84")
  marked_points[, 1:2] <- coordinates(spTransform(xy, CRS("+proj=utm +zone=53 ellps=WGS84")))
  
  marked_points$tp <- (marked_points$latitude - min(marked_points$latitude))/(max(marked_points$latitude) - min(marked_points$latitude))
  owin <- owin(xrange = c(mean(marked_points$longitude), mean(marked_points$longitude)+500), yrange = c(min(marked_points$latitude), max(marked_points$latitude)))
  transect_bounds <- ppp(x = rep(mean(marked_points$longitude), 47*2), y = rep(c(min(marked_points$latitude), max(marked_points$latitude)), 47), window = owin)
  jointed_vertices <- matrix(FALSE, nrow = 94, ncol = 94)
  pairs <- cbind(seq(1, 93), seq(2, 94))
  jointed_vertices[rbind(pairs, cbind(pairs[, 2], pairs[, 1]))] <- TRUE
  
  
  transect_bounds <- ppp(x = rep(mean(marked_points$longitude), 2), y = rep(c(min(marked_points$latitude), max(marked_points$latitude)), 1), window = owin)
  transect_line <- linnet(vertices = transect_bounds, m = jointed_vertices[1:2, 1:2])
  marked_points$longitude <- mean(transect_bounds$x)
  names(marked_points) <- c("x", "y", "mark", "seg", "tp")
  points_day <- marked_points[marked_points$seg == day, ]
  points_day$seg <- 1
  point_network <- lpp(X = points_day, L = transect_line)
  pair_correlation <- linearpcfcross(point_network, "fish", "dolphin")
  
  r <- c(r, pair_correlation$r)
  est <- c(est, pair_correlation$est)

}

cor_mean <- aggregate(est, list(r), mean)
plot(cor_mean$Group.1, cor_mean$x)



#linear K cross
K <- linearKcross(point_network, "net", "dolphin")
plot(K)

#pair correlation function
pair_correlation <- linearpcfcross(point_network, "net", "dolphin")
plot(pair_correlation) 


#adding sandy beaches without shark nets
K <- linearKcross(point_network, "no_net", "dolphin")
plot(K)

pair_correlation <- linearpcfcross(point_network, "no_net", "dolphin")
plot(pair_correlation) 

#combining all beaches together
marked_points$mark <- revalue(marked_points$mark, c(no_net = "beach", net = "beach"))
point_network <- lpp(X = marked_points, L = transect_line) #check duplicated values

pair_correlation <- linearpcfcross(point_network, "beach", "dolphin")
plot(pair_correlation) 


#is a sighting closer to a netted beach or non-netted beach


#closest net
closest_net   <- apply(sighting_to_net, 1, min, na.rm = TRUE)

#closest no_net
dist_to_no_net  <- distAllCombinations(dat_bot$Lat, no_nets$latitude, no_nets$longitude)
closest_no_net   <- apply(dist_to_no_net, 1, min, na.rm = TRUE)

closest_to_netted <- closest_net < closest_no_net
closest_to_netted[closest_net > 1000 & closest_no_net > 1000] <- NA #can specify maximum distance from net/no_net


table(closest_to_netted) #78% of sightings are closer to a netted beach than a non-netted beach, however 77% of beaches are netted

#77% of beaches are netted
#78% of sightings are closer to a netted beach than a non-netted beach

#84% of sightings withing 1km of a beach are at a netted beach vs non-netted


