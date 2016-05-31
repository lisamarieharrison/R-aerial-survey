#aerial survey summary for report
#author: Lisa-Marie Harrison
#date: 22/04/2016

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
library(splines)
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


#remove duplicated observations
#dat <- dat[!duplicated(dat), ]

#remove secondary observations
dat <- dat[!dat$Secondary == "Y", ]

#check number of sightings = 2695
#getting 2375 
dim(dat[dat$Type == "S", ])

dat$Count <- NA
dat$Count[dat$Type == "S"] <- 1

dat$Date <- chron(dates. = as.character(dat$Date), format = "d/m/y")

#plot dolphin group size
hist(dat$Number[dat$Species == "BOT"], main = "Bottlenose Dolphin Group Size", xlab = "Number of dolphins", col = "lightgrey")

#combine all non-hammerhead sharks
dat$Species[dat$Species %in% c("W", "Wh", "BS")] <- "S"


seasonFromDate <- function(x) {
  
  season_lookup <- data.frame("month" = 1:12, "season" = c(rep("Summer", 2), rep("Autumn", 3), rep("Winter", 3), rep("Spring", 3), "Summer"))
  season <- season_lookup[as.numeric(months(chron(dates. = x, format = "d/m/y"))), "season"]
  season <- factor(season)
  
  return (season)
  
}


#sightings/survey by season
sightings_per_survey <- table(dat$Species, dat$Date)
season <- seasonFromDate(colnames(sightings_per_survey))

layout.show(layout(matrix(c(1,2,0,0,3,4),3,2,byrow=TRUE), widths=c(1,1), heights=c(4,0.01,4)))
par(bty = "n", oma = c(0, 4, 0, 0))
boxplot(sightings_per_survey["B", ] ~ season, col = c("indianred", "skyblue", "palegoldenrod"), main = "Baitfish")
boxplot(sightings_per_survey["BOT", ] ~ season, col = c("indianred", "skyblue", "palegoldenrod"), main = "Bottlenose dolphin")
boxplot(sightings_per_survey["S", ] ~ season, col = c("indianred", "skyblue", "palegoldenrod"), main = "Sharks")
boxplot(sightings_per_survey["HH", ] ~ season, col = c("indianred", "skyblue", "palegoldenrod"), main = "Hammerhead Sharks")
mtext('Sightings per survey', side = 2, outer = TRUE, line = 2)


#baitfish size by distance (south only)
x <- revalue(dat$Length.Size, c("S"="small", "M"="med", "L" = "large"))
par(mfrow = c(1, 2), bty = "n", oma = c(1, 3, 1, 0), mar = c(5, 1, 5, 2))
boxplot(dat$Dist..from.transect[dat$Species == "B" & dat$Length.Size !="" & dat$Flight.Direction == "N" & dat$Dist..from.transect <= 300] ~ factor(x[dat$Species == "B" & dat$Length.Size !="" & dat$Flight.Direction == "N"& dat$Dist..from.transect <= 300], c("small", "med", "large")), 
        pch = 19, col = c("indianred", "skyblue", "palegoldenrod"), main = "North")
boxplot(dat$Dist..from.transect[dat$Species == "B" & dat$Length.Size !="" & dat$Flight.Direction == "S"] ~ factor(x[dat$Species == "B" & dat$Length.Size !="" & dat$Flight.Direction == "S"], c("small", "med", "large")), 
        pch = 19, col = c("indianred", "skyblue", "palegoldenrod"), main = "South")
mtext("Distance from transect (m)", side = 2, outer = TRUE, line = 2)


#distance offshore
par(mfrow = c(1, 1), mar = c(6, 5, 1, 1), oma = c(0, 0, 0, 0))
dat_south_boxplot <- dat[dat$Flight.Direction == "S" & dat$Species %in% c("B", "BOT", "S", "HH", "HB", "R", "T", "P"), ]
dat_south_boxplot$Species <- factor(dat_south_boxplot$Species)
boxplot(dat_south_boxplot$Dist..from.transect ~ dat_south_boxplot$Species, ylab = "Distance (m)", ylim = c(50, 1000), xaxt = "n", col = "lightgrey", pch = 19)
text(1:8, y=par()$usr[3]-0.1*100,
     labels=c("Baitfish", "Bottlenose Dolphin", "Humpback Whale", "Hammerhead Shark", "Seal", "Ray", "Shark", "Turtle"), srt=45, adj=1, xpd=TRUE)

# ----------------------------------- ENVIRONMENTAL FACTORS -------------------------------------- #

dat$Time <- chron(times. = dat$Time, format = "h:m:s")
dat <- dat[!is.na(dat$Time), ]

dat_north <- dat[dat$Flight.Direction == "N", ]
dat_south <- dat[dat$Flight.Direction == "S", ]

#calculate number of hours flown in each direction
effort_north <- calculateEffort(dat_north)
envt_effort <- effort_north$envt.var
hours_flown_north <- effort_north$total/60


effort_south <- calculateEffort(dat_south)
envt_effort_south <- effort_south$envt.var
hours_flown_south <- effort_south$total/60

dat_season_north <- dat_north[dat_north$Species %in% c("B", "BOT", "HH", "S"), ]
dat_season_south <- dat_south[dat_south$Species %in% c("B", "BOT", "HH", "S"), ]
dat_season_north$Species <- factor(dat_season_north$Species)
dat_season_south$Species <- factor(dat_season_south$Species)

#proportion of time and each environmental factor level

par(mfrow = c(2, 2), mar = c(2, 2, 2, 3))
sea_state_tab <- table(dat$Beaufort.Sea.State, dat$Season)
season_totals <- colSums(sea_state_tab)
weighted_sea_state <- sea_state_tab/rep(season_totals, each = nrow(sea_state_tab))
barplot(weighted_sea_state, legend = T, main = "Sea State", args.legend = list(x = ncol(weighted_sea_state) + 1.5, y=1, bty = "n", cex = 0.8))

turbidity_tab <- table(dat$Water.clarity, dat$Season)
season_totals <- colSums(turbidity_tab)
weighted_turbidity <- turbidity_tab/rep(season_totals, each = nrow(turbidity_tab))
barplot(weighted_turbidity, legend = T, main = "Turbidity", args.legend = list(x = ncol(weighted_sea_state) + 1.5, y=1, bty = "n", cex = 0.8)) 

cloud_tab <- table(dat$Cloud.cover, dat$Season)
season_totals <- colSums(cloud_tab)
weighted_cloud <- cloud_tab/rep(season_totals, each = nrow(cloud_tab))
barplot(weighted_cloud, legend = T, main = "Cloud Cover", args.legend = list(x = ncol(weighted_sea_state) + 1.5, y=1, bty = "n", cex = 0.8)) 


#sightings at each environmental variable level


par(mfrow = c(1, 2))

#sea state
north_tab <- table(dat_season_north$Beaufort.Sea.State, dat_season_north$Species)
species <- which(colnames(north_tab) %in% c("B", "BOT", "HH", "W"))
north_tab_effort <- north_tab[, species]/(envt_effort[c(1, 4, 6)]/60)

barplot(north_tab_effort, legend = T, main = "Sightings/hr at each sea state", args.legend = list(x = 5, y=20, bty = "n"))

#turbidity
north_tab <- table(dat_season_north$Water.clarity, dat_season_north$Species)
north_tab_effort <- north_tab[, species]/(envt_effort[c(17:19)]/60)

barplot(north_tab_effort, main = "Sightings/hr at each turbidity", legend = T, args.legend = list(x = 5, y=20, bty = "n"))


par(mfrow = c(1, 2))

#sea state
south_tab <- table(dat_season_south$Beaufort.Sea.State, dat_season_south$Species)
species <- which(colnames(south_tab) %in% c("B", "BOT", "HH", "S"))
south_tab_effort <- south_tab[, species]/(envt_effort_south[, paste0("Beaufort.Sea.State.", sort(unique(dat_season_south$Beaufort.Sea.State)))]/60) 

barplot(south_tab_effort, legend = T, main = "Sea State", args.legend = list("topright", bty = "n"), cex.names = 0.8)

#turbidity
south_tab <- table(dat_season_south$Water.clarity, dat_season_south$Species)
south_tab_effort <- south_tab/(envt_effort_south[, paste0("Water.clarity.", sort(unique(dat_season_south$Water.clarity)))]/60) 

barplot(south_tab_effort, main = "Turbidity", legend = T, args.legend = list("topright", bty = "n"), cex.names = 0.8)

# --------------------------------------- DISTANCE SAMPLING --------------------------------------------#


lisa_obs <- dat[dat$Observer == "Lisa" & dat$Flight.Direction == "S", ]

#add a unique trial number to each day
lisa_obs$Trial <- as.numeric(as.factor(as.character(lisa_obs$Date)))


#plots of detection probability for each species

par(mfrow = c(2, 2))

#Baitfish
total_observations <- createDistanceData(species = "B", lisa_obs, truncate = 1000, direction = "S")
det_fun_B <- ddf(method = 'ds',dsmodel =~ cds(key = "gamma", formula=~1), data = total_observations, meta.data = list(left = 50, width = 500))
plot(det_fun_B, main = "Baitfish")

#Bottlenose dolphins
total_observations <- createDistanceData(species = "BOT", lisa_obs, truncate = 1000, direction = "S")
det_fun_BOT <- ddf(method = 'ds',dsmodel =~ cds(key = "gamma", formula=~1), data = total_observations, meta.data = list(left = 50, width = 500))
plot(det_fun_BOT, main = "Bottlenose Dolphins")

#Hammerhead
total_observations <- createDistanceData(species = "HH", lisa_obs, truncate = 1000, direction = "S")
det_fun_HH <- ddf(method = 'ds',dsmodel =~ cds(key = "hr", formula=~1), data = total_observations, meta.data = list(left = 50, width = 1000))
plot(det_fun_HH, main = "Hammerhead Sharks")

#Other Sharks
total_observations <- createDistanceData(species = "S", lisa_obs, truncate = 1000, direction = "S")
det_fun_S <- ddf(method = 'ds',dsmodel =~ cds(key = "hr", formula=~1), data = total_observations, meta.data = list(left = 50, width = 1000))
plot(det_fun_S, main = "Sharks")


#abundance calculations
#using p(0) estimated from interobserver surveys

calcAbundanceAndCV(det_fun_B, line_length = 265, n_surveys = 47, p_0 = 0.647, p_0_cv = 0.025)

calcAbundanceAndCV(det_fun_BOT, line_length = 265, n_surveys = 47, p_0 = 0.778, p_0_cv = 0.057)
calcAbundanceAndCV(det_fun_BOT, line_length = 265, n_surveys = 47, p_0 = 0.778, p_0_cv = 0.057, group = FALSE)

calcAbundanceAndCV(det_fun_S, line_length = 265, n_surveys = 47, p_0 = 0.500, p_0_cv = 0.500)


#------------------------------------- INTER-OBSERVER -------------------------------------------------#


dat <- read.csv("interobserver_20150522.csv", header = T)


lisa <- dat[dat$Observer == "Lisa", ]
vic <- dat[dat$Observer == "Vic", ]


checkSame <- function(r1, r2) {
  #checks if 2 observations are the same and returns boolean
  if (r1$Date == r2$Date & substr(r1$Time, 1, 5) == substr(r2$Time, 1, 5) & r1$Species == r2$Species) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}



findOverlap <- function(set1, set2, i, j,overlap, lisa_missed, vic_missed) {
  
  
  if (j > nrow(set2)) {
    #print(paste(i, "No overlapping observation found. Try next observation"))
    lisa_missed <- rbind(lisa_missed, set1[i, ])
    return(list(overlap = overlap, vic_missed = vic_missed, lisa_missed = lisa_missed))
  }
  
  if (checkSame(set1[i, ], set2[j, ])) {
    overlap <- rbind(overlap, set1[i, ])
    #print("Success: found a same observation!")
    vic_missed <- vic_missed[-which(vic_missed$Date == set1$Date[i] & substr(vic_missed$Time, 1, 5) == substr(set1$Time[i], 1, 5) & vic_missed$Species == set1$Species[i])[1], ]
    return(list(overlap = overlap, vic_missed = vic_missed, lisa_missed = lisa_missed))
  } else {
    findOverlap(set1, set2, i = i, j = j + 1, overlap = overlap, lisa_missed = lisa_missed, 
                vic_missed = vic_missed)
  }
  
}

overlap <- dat[0, ]
vic_missed <- lisa
lisa_missed <- dat[0, ]

for (i in 1:nrow(vic)) {
  set2 <- vic_missed[vic_missed$Date == vic$Date[i], ]
  
  result <- findOverlap(set1 = vic, set2 = set2, i = i, j = 1, 
                        overlap = overlap, lisa_missed = lisa_missed, 
                        vic_missed = vic_missed)
  overlap <- result$overlap
  lisa_missed <- result$lisa_missed
  vic_missed <- result$vic_missed
}
lisa_missed <- lisa_missed[-1, ]

total_missed <- rbind(lisa_missed, vic_missed)
total_missed$Species[total_missed$Species %in% c("BS", "W", "Wh")] <- "S"

#missed sightings per survey

total_interobs <- rbind(overlap, lisa_missed, vic_missed)
total_B <- colSums(table(total_interobs$Observer[total_interobs$Species == "B"], total_interobs$Date[total_interobs$Species == "B"]))

par(mfrow = c(1, 2))
hist(table(total_missed$Observer[total_missed$Species == "B"], total_missed$Date[total_missed$Species == "B"])[1, ]/total_B, 
     main = "Vic", xlab = "", col = "lightgrey", xlim = c(0, 1))
hist(table(total_missed$Observer[total_missed$Species == "B"], total_missed$Date[total_missed$Species == "B"])[2, ]/total_B, 
     main = "Lisa", xlab = "", col = "lightgrey", xlim = c(0, 1))


total_obs <- rbind(lisa_missed, vic_missed, overlap, lisa_missed, vic_missed, overlap)
total_obs$observer <- c(rep(1, nrow(total_obs)/2), rep(2, nrow(total_obs)/2))
total_obs$detected <- c(rep(0, nrow(lisa_missed)), rep(1, nrow(vic_missed)), rep(1, nrow(overlap)),
                        rep(1, nrow(lisa_missed)), rep(0, nrow(vic_missed)), rep(1, nrow(overlap)))
total_obs$object <- rep(1:(nrow(total_obs)/2), 2)
total_obs$distance <- total_obs$Dist..from.transect
total_obs <- total_obs[!is.na(total_obs$distance), ]


missed_table <- table(total_obs$observer[total_obs$detected == 0], total_obs$Species[total_obs$detected == 0])
missed_table <- missed_table[, c("B", "BOT", "HH", "S", "R", "T")]
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(8, 5, 1, 1))
barplot(missed_table, beside = TRUE, legend = c("Observer 1", "Observer 2"), ylab = "Number of missed sightings", args.legend = list(bty = "n"), xaxt = "n")
text(c(2, 5, 8, 11, 14, 17), y=par()$usr[3]-0.1*100,
     labels=c("Baitfish", "Bottlenose Dolphin", "Hammerhead Shark", "Shark", "Ray", "Turtle"), adj=1, xpd=TRUE, srt = 45)


#add a unique trial number to each day
total_obs$Trial <- as.numeric(as.factor(as.character(total_obs$Date)))
total_obs$count <- 1

interobs_aggregate <- aggregate(count ~ observer + Species + Trial,data = total_obs[total_obs$detected == 0, ], FUN=sum)
interobs_aggregate$Species[interobs_aggregate$Species %in% c("W", "BS", "Wh")] <- "S"
interobs_aggregate <- interobs_aggregate[interobs_aggregate$Species %in% c("B", "BOT", "S"), ]
interobs_aggregate$Species <- factor(interobs_aggregate$Species)
total_obs$Species[total_obs$Species %in% c("W", "BS", "Wh")] <- "S"
sightings_per_trial <- table(total_obs$Species[!duplicated(total_obs$object)], total_obs$Trial[!duplicated(total_obs$object)])
interobs_aggregate$total_observations <- diag(sightings_per_trial[match(interobs_aggregate$Species, rownames(sightings_per_trial)), match(interobs_aggregate$Trial, colnames(sightings_per_trial))])

#add no missed sightings in
for (species in c("B", "BOT", "S")) {
  
    for (observer in 1:2) {
      
      missing <- setdiff(1:18, as.numeric(names(table(interobs_aggregate$Trial[interobs_aggregate$Species == species & interobs_aggregate$observer == observer]))))
      
      if (length(missing) > 0) {
        
        interobs_aggregate <- rbind(interobs_aggregate, cbind("observer" = rep(observer, length(missing)), "Species" = rep(species, length(missing)), "Trial" =  missing, 
                                                              "count" = rep(0, length(missing)), "total_observations" =  sightings_per_trial[species, missing]))
        
      }
      
    }
  
}
interobs_aggregate$count <- as.numeric(interobs_aggregate$count)
interobs_aggregate$total_observations <- as.numeric(interobs_aggregate$total_observations)
interobs_aggregate$Trial <- as.numeric(interobs_aggregate$Trial)

boxplot(count/total_observations ~ observer + Species, data = interobs_aggregate, ylab = "Proportion of missed sightings per survey", col = c("darkgrey", "lightgrey"), pch = 19, xaxt = "n")
axis(1, at = 1:6, rep(c(1, 2), 3))
axis(1, at = c(1.5, 3.5, 5.5), c("Baitfish", "Bottlenose dolphins", "Sharks"), pos = -0.1, line = F, tick = F)
axis(1, at = 3.5, c("Observer and Species"), pos = -0.2, line = F, tick = F)


#learning curve

interobs_aggregate <- interobs_aggregate[order(interobs_aggregate$Trial), ]

species_lookup <- data.frame("Species" = c("B", "BOT", "S"), "Full_name" = c("Baitfish", "Bottlenose dolphins", "Sharks"))

par(mfrow = c(1, 2), oma = c(4, 4, 6, 0), mar = c(1, 1, 2, 1))

for (s in c("B", "BOT")) {
  
  plot(interobs_aggregate$Trial[interobs_aggregate$observer == 1 & interobs_aggregate$Species == s], interobs_aggregate$count[interobs_aggregate$observer == 1 & interobs_aggregate$Species == s]/interobs_aggregate$total_observations[interobs_aggregate$observer == 1 & interobs_aggregate$Species == s], type = "l", xlab = "", ylab = "", ylim = c(0, 1), main = species_lookup$Full_name[species_lookup$Species == s])
  points(interobs_aggregate$Trial[interobs_aggregate$observer == 2 & interobs_aggregate$Species == s], interobs_aggregate$count[interobs_aggregate$observer == 2 & interobs_aggregate$Species == s]/interobs_aggregate$total_observations[interobs_aggregate$observer == 2 & interobs_aggregate$Species == s], type = "l", col = "red")
  mtext('Trial', side = 1, outer = TRUE, line = 2)
  mtext('Proportion of missed sightings', side = 2, outer = TRUE, line = 2)
  
}
legend("topright", c("observer 1", "observer 2"), col = c("black", "red"), lwd = 2, bty = "n", cex = 0.75)


# distance sampling with gamma detection function


region.table <- data.frame("Region.Label" = 1, "Area" = 265)
sample.table <- data.frame("Sample.Label" = 1:46, "Region.Label" = 1, "Effort" = 265)
obs.table    <- data.frame("object" = 1:nrow(total_observations), "Region.Label" = 1, "Sample.Label" = total_observations$Trial)

det_fun_S <- ds(data = total_observations, truncation = list(left = 50, right = 1000), key = "hr", region.table=region.table, sample.table=sample.table, obs.table=obs.table, convert.units=0.001)

distance.sample.size(cv.pct = 30, N = 100, detection.function = "hazard", theta = c(1.640992, 5.508635), w = 1)


p_total <- ddf(method = 'trial',dsmodel =~ cds(key = "gamma", formula=~1), mrmodel =~ glm(formula=~1),
               data = total_obs[total_obs$Species == "BOT", ], meta.data = list(left = 50, width = 500))


p_total <- ddf(method = 'trial',dsmodel =~ cds(key = "gamma", formula=~1), mrmodel =~ glm(formula=~bs(distance,degree=4)),
               data = total_obs[total_obs$Species == "BOT", ], meta.data = list(left = 50, width = 500))

#find apex of ds model
apex <- mrds:::apex.gamma(det_fun_BOT$ds$aux$ddfobj)[1] + 50 #need to add 50 because data left truncated

plot(p_total$mr)
abline(v = apex + 50, col = "red")


#from summary.trial.fi
#substitute apex for 0 to get the p(0) value returned in summary(ddf.obj)

newdat <- p_total$data
newdat <- newdat[newdat$distance <= p_total$meta.data$width &
                   newdat$distance >= p_total$meta.data$left, ]
newdat <- newdat[newdat$observer == 1 & newdat$detected == 1, ]

newdat$distance <- rep(apex,length(newdat$distance))


pred.at0 <- predict(p_total$mr,newdat,type="response")[[1]]
pdot <- p_total$mr$fitted[, 1]
avgp <- function(model,pdot,...) {return(pdot)}
n <- length(newdat$distance)/2
vcov <- solvecov(p_total$mr$hessian)$inv
average.p <-  length(newdat$distance[newdat$observer == 1 & newdat$detected == 1])/p_total$mr$Nhat
average.p0.1 <- sum(avgp(p_total$mr,pred.at0)/pdot)
NCovered = NCovered.trial(p_total$mr$par, p_total$mr)
var.pbar.list <- prob.se(p_total$mr,avgp,vcov,fittedmodel=NULL)
se.obj <- calc.se.Np(p_total$mr, avgp, n, average.p)
Nhatvar.list <- se.obj$Nhatvar.list
cvN <- se.obj$Nhat.se/p_total$mr$Nhat
covar <- t(Nhatvar.list$partial) %*% vcov %*% var.pbar.list$partial +
  var.pbar.list$covar
var.pbar <- (average.p0.1/p_total$mr$Nhat)^2*(cvN^2 + var.pbar.list$var/average.p0.1^2
                                              - 2*covar/(average.p0.1*p_total$mr$Nhat))

#p(y) where y = apex and corresponding se
data.frame("average.p0.1" = average.p0.1/p_total$mr$Nhat, "se" = sqrt(var.pbar))



