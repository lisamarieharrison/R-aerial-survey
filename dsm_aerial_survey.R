#dsm for comparison with mrds and mcmc aerial survey analysis
#author: Lisa-Marie Harrison
#date: 6/9/2016

dat <- read.csv("C:/Users/43439535/Documents/Lisa/phd/aerial survey/data/aerial_survey_summary_r.csv", header = T)

if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/aerial survey/data")
  source_location <- "~/Lisa/phd/aerial survey/R/R-aerial-survey/"
} else {
  setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/data")
  source_location <- "~/phd/aerial survey/R/code/R-aerial-survey/"
}

library(chron)
library(mrds)
library(dsm)
library(sp)
library(rgeos) #gLength

dat <- read.csv("aerial_survey_summary_r.csv", header = T)

source_list <- c("functions/createDistanceData.R",
               "functions/gamma_det_fun.R")

invisible(Map(source, paste0(source_location, source_list)))
source("C:/Users/43439535/Documents/Lisa/phd/Mixed models/R code/R-functions-southern-ocean/grid_plot_obj.R")

#remove secondary observations
dat <- dat[!dat$Secondary == "Y", ]

#check number of sightings = 2695
#getting 2375 
dim(dat[dat$Type == "S", ])
dat$Count <- NA
dat$Count[dat$Type == "S"] <- 1
dat$Date <- chron(dates. = as.character(dat$Date), format = "d/m/y")

#combine all non-hammerhead sharks
dat$Species[dat$Species %in% c("W", "Wh", "BS")] <- "S"

lisa_obs <- dat[dat$Observer == "Lisa" & dat$Flight.Direction == "S", ]

#add a unique trial number to each day
lisa_obs$Trial <- as.numeric(as.factor(as.character(lisa_obs$Date)))

#Fit detection function
#same as for mrds because dsm ds doesn't support gamma detection function
total_observations <- createDistanceData(species = "BOT", lisa_obs, truncate = 1000, direction = "S")
det_fun <- ddf(method = 'ds',dsmodel =~ mcds(key = "gamma", formula=~1), data = total_observations, meta.data = list(left = 0, width = 1000))
plot(det_fun)

#transform observation coordinates to northings
total_observations_latlon <- SpatialPoints(apply(total_observations[, c("long", "lat")], 2, as.numeric), proj4string = CRS("+proj=longlat +datum=WGS84"))
total_observations_utm <- spTransform(total_observations_latlon, CRS("+proj=utm +zone=53 ellps=WGS84"))
total_observations$E <- coordinates(total_observations_utm)[, 1]
total_observations$N <- coordinates(total_observations_utm)[, 2]
total_observations$N[total_observations$N > 0] <- -1*total_observations$N[total_observations$N > 0]

# ------------------------------ CREATE DATA FRAMES FOR DSM ------------------------------------- #


#read cruise track
track <- readGPX("C:/Users/43439535/Documents/Lisa/phd/aerial survey/data/gps data/tracks/cruise_track.gpx")
track_coords <- track$tracks[[2]]$`ACTIVE LOG`

#make spatial Lines object
track_line <- SpatialLines(list(Lines(list(Line(cbind(track_coords$lon, track_coords$lat))), ID = 1)), proj4string = CRS("+proj=longlat +datum=WGS84"))

#transform to utm
track_utm <- spTransform(track_line, CRS("+proj=utm +zone=53 ellps=WGS84"))

#find total length of track
track_length <- gLength(track_utm)

#break track into regularly sized segments
segment_size <- 4000 #m
track_segmented <- spsample(track_utm, track_length/segment_size, type = "regular")

#find segment mid points
mid_points <- coordinates(track_segmented)
transects <- sort(unique(total_observations$Trial))

#create segdata for dsm
seg_data <- data.frame("Effort" = segment_size, "Seg.Label" = rep(1:nrow(mid_points), length(transects)), "Transect.Label" = rep(transects, each = nrow(mid_points)),
                       "X" = mid_points[, 1], "Y" = mid_points[, 2])


#add baitfish to correct segments
seg_data$fish <- 0

fish_observations <- createDistanceData(species = "B", lisa_obs, truncate = 1000, direction = "S")
fish_observations$segment <- NA
fish_observations_latlon <- SpatialPoints(apply(fish_observations[, c("long", "lat")], 2, as.numeric), proj4string = CRS("+proj=longlat +datum=WGS84"))
fish_observations_utm <- spTransform(fish_observations_latlon, CRS("+proj=utm +zone=53 ellps=WGS84"))
fish_observations$E <- coordinates(fish_observations_utm)[, 1]
fish_observations$N <- coordinates(fish_observations_utm)[, 2]
fish_observations$N[fish_observations$N > 0] <- -1*fish_observations$N[fish_observations$N > 0]

for (obs in 1:nrow(fish_observations)) {
  
  dists <- abs(fish_observations$E[obs] - mid_points[, 1]) + abs(fish_observations$N[obs] - mid_points[, 2])
  fish_observations$segment[obs] <-  which.min(dists)
  
}

fish_observations$Sample.Label <- paste0(fish_observations$Trial, "-", fish_observations$segment)
seg_data$fish[na.omit(match(names(table(fish_observations$Sample.Label)), seg_data$Sample.Label))] <- table(fish_observations$Sample.Label)[!is.na(match(names(table(fish_observations$Sample.Label)), seg_data$Sample.Label))]


#find which segment each observation is closest to
total_observations$segment <- NA
total_observations$dist_to_seg <- NA #distance to closest segment for error checking
for (i in 1:nrow(total_observations)) {
  dists <- abs(total_observations$E[i] - mid_points[, 1]) + abs(total_observations$N[i] - mid_points[, 2])
  total_observations$segment[i] <-  which.min(dists)
  total_observations$dist_to_seg[i] <- min(dists)
}

#create obsdata for dsm
obs_data <- data.frame("object" = total_observations$object, "Seg.Label" = total_observations$segment, "size" = total_observations$size, 
                       "distance" = total_observations$distance, "Transect.Label" = total_observations$Trial)

obs_data$Sample.Label <- paste0(obs_data$Transect.Label, "-", obs_data$Seg.Label)
seg_data$Sample.Label <- paste0(seg_data$Transect.Label, "-", seg_data$Seg.Label)



bot_dsm <- dsm(count ~ s(X, Y) + s(fish), ddf.obj = det_fun, segment.data = seg_data, observation.data = obs_data)

summary(bot_dsm)
plot(bot_dsm)


#predict count
#uses mean number of fish schools per segment
pred_dat <- data.frame(cbind(mid_points, aggregate(seg_data$fish, list(seg_data$Seg.Label), mean)$x))
colnames(pred_dat) <- c("X", "Y", "fish")
pred_num <- predict(bot_dsm, newdata = pred_dat, off.set = segment_size*1000)


qplot(pred_dat[, 1], pred_dat[, 2], colour = pred_num, size = 4, xlab = "Easting", ylab = "Northing")

#total estimate of groups
sum(pred_num)

