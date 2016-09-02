
dat <- read.csv("C:/Users/43439535/Documents/Lisa/phd/aerial survey/data/aerial_survey_summary_r.csv", header = T)

if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/aerial survey/data")
  source_location <- "~/Lisa/phd/aerial survey/R/R-aerial-survey/"
} else {
  setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/data")
  source_location <- "~/phd/aerial survey/R/code/R-aerial-survey/"
}

library(chron)
library(plyr)
library(mrds)
library(Distance)

dat <- read.csv("aerial_survey_summary_r.csv", header = T)

file_list <- c("fun_calculate_envt_effort.R",
               "functions/createDistanceData.R",
               "/functions/mrds_modified_functions.R",
               "functions/calcAbundanceAndCVExtraP.R",
               "mrds_modified_for_gamma.R")

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
total_observations <- createDistanceData(species = "B", lisa_obs, truncate = 1000, direction = "S")
det_fun <- ddf(method = 'ds',dsmodel =~ mcds(key = "gamma", formula=~1), data = total_observations, meta.data = list(left = 0, width = 1000))
plot(det_fun)

est <- exp(det_fun$par) + c(1, 0)
sd <- c(summary(det_fun)$coeff$key.shape[2]$se, summary(det_fun)$coeff$key.scale[2]$se)
cv <- sd/est

#summary table

sum_sats <- cbind(est, sd, cv)
colnames(sum_sats) <- c("Est", "SE", "CV")
rownames(sum_sats) <- c("Shape", "Scale")
sum_sats <- round(sum_sats, 4)
if (det_fun$ds$aux$ddfobj$type == "hn") { #if half normal model, only scale parameter estimated
  sum_sats <- sum_sats[2, ]
}

sum_sats

#average p
summary(det_fun)$average.p


#calculate abundance
calcAbundanceAndCV(det_fun, line_length = 265, n_surveys = length(unique(total_observations$Trial)))

a <- dht(det_fun, region.table = data.frame("Region.Label" = 1, "Area" = 265000), sample.table = data.frame("Region.Label" = 1, "Sample.Label" = 1:46, "Effort" = 265),
      obs.table = data.frame("object" = 1:nrow(total_observations), "Region.Label" = 1, "Sample.Label" = total_observations$Trial))


