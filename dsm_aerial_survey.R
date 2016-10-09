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
library(raster)
library(maps)
library(mapdata)
library(maptools)

dat <- read.csv("aerial_survey_summary_r.csv", header = T)

file_list <- c("functions/createDistanceData.R",
               "gamma_det_fun.R")

for (f in file_list) {
  
  source(paste0(source_location, f))
  
}

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

obs_data <- data.frame("object" = total_observations$object, "Sample.Label" = total_observations$Trial, "size" = total_observations$size, 
                       "distance" = total_observations$distance, "Transect.Label" = total_observations$Trial)

seg_data <- data.frame("Effort", "Sample.Label", "Transect.Label" = sort(unique(total_observations$Trial)))


#map of coastline


lat_bound <- c(-34.6, -32.6)
long_bound <- c(150.8, 152)

aust  <- map("world2Hires", regions="Australia", fill = TRUE, col = "grey", xlim = long_bound, ylim = lat_bound, plot = FALSE)
aust_poly <- map2SpatialPolygons(aust, IDs = aust$names, proj4string=CRS("+proj=longlat +datum=WGS84"))

syd_poly <- crop(aust_poly, extent(c(long_bound, lat_bound)))

#make empty raster
grid <- raster(extent(aust_poly))

# Choose its res (m)
res(grid) <- 1

plot(grid)

