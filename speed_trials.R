dat <- read.csv("C:/Users/43439535/Documents/Lisa/phd/aerial survey/data/speed trials/completed datasheets/Completed_matched_data.csv")
lisa <- read.csv("C:/Users/43439535/Documents/Lisa/phd/aerial survey/data/speed trials/completed datasheets/speed_trials_datasheet_lisa.csv", skip = 1)
sally <- read.csv("C:/Users/43439535/Documents/Lisa/phd/aerial survey/data/speed trials/completed datasheets/speed_trials_datasheet_Sally.csv", skip = 1)
vic <- read.csv("C:/Users/43439535/Documents/Lisa/phd/aerial survey/data/speed trials/completed datasheets/speed_trials_datasheet_vic.csv", skip = 1)
rob <- read.csv("C:/Users/43439535/Documents/Lisa/phd/aerial survey/data/speed trials/completed datasheets/speed_trials_datasheet_Rob.csv", skip = 1)

h1_speed <- c(80, 100, 100, 80, 80, 100, 100, 80, 80, 100, 100, 80, 80, 100, 100, 80)
h2_speed <- c(100, 80, 80, 100, 100, 80, 80, 100, 100, 80, 80, 100, 100, 80, 80, 100)

dat[dat == 0] <- NA

dat$distance <- rowMeans(cbind(lisa$Dist_from_transect[dat$Lisa - 2], sally$Dist_from_transect[dat$Sally - 2], vic$Dist_from_transect[dat$Vic - 2], rob$Dist_from_transect[dat$Rob - 2]), na.rm = TRUE)


dat<- rbind(dat, dat, dat, dat)
dat$observer <- c(rep(1, nrow(dat)/4), rep(2, nrow(dat)/4), rep(3, nrow(dat)/4), rep(4, nrow(dat)/4))
dat$detected <- 0
dat$object <- dat$Unique_objects

observer_lookup <- data.frame("observer" = c("Lisa", "Sally", "Vic", "Rob"), "number" = 1:4)

for (i in 1:4) {
  
  obs <- dat[, as.character(observer_lookup[i, 1])]
  
  dat$detected[dat$observer == i & !is.na(dat[, as.character(observer_lookup[i, 1])])] <- 1

}

dat$speed[dat$observer %in% 1:2] <- h1_speed[dat$Transect][dat$observer %in% 1:2]
dat$speed[dat$observer %in% 3:4] <- h2_speed[dat$Transect][dat$observer %in% 3:4]

dat$direction <- "N"
dat$direction[dat$Transect %% 2 == 0] <- "S"

dat <- dat[!dat$object %in% which(table(dat$object) == 8), ]
dat <- dat[!is.nan(dat$distance), ]

h1 <- dat[dat$observer %in% 1:2 & !(is.na(dat$Sally) & is.na(dat$Lisa)), ]
h2 <- dat[dat$observer %in% 3:4 & !(is.na(dat$Vic) & is.na(dat$Rob)), ]

combined <- rbind(h1, h2)
combined$team <- 1
combined$team[(nrow(h1) + 1):nrow(combined)] <- 2

b_80 <- combined[combined$Species == "B" & combined$speed == 80, ]
b_80$object[!(is.na(b_80$Sally) & is.na(b_80$Lisa)) & b_80$observer %in% 1:2] <- rep(1:(nrow(b_80[!(is.na(b_80$Sally) & is.na(b_80$Lisa)) & b_80$observer %in% 1:2, ])/2), 2)
b_80$object[!(is.na(b_80$Vic) & is.na(b_80$Rob)) & b_80$observer %in% 3:4] <- rep((nrow(b_80[!(is.na(b_80$Sally) & is.na(b_80$Lisa)) & b_80$observer %in% 1:2, ])/2 + 1):((nrow(b_80[!(is.na(b_80$Sally) & is.na(b_80$Lisa)) & b_80$observer %in% 1:2, ])/2) + (nrow(b_80[!(is.na(b_80$Vic) & is.na(b_80$Rob)) & b_80$observer %in% 3:4, ])/2)), 2)

b_80$observer[b_80$observer == 3] <- 1
b_80$observer[b_80$observer == 4] <- 2

b_80$Region.Label <- 1
b_80$Sample.Label <- b_80$Transect


det_b_80 <- ddf(method = 'io',dsmodel =~ cds(key = "gamma", formula=~team + direction), mrmodel =~ glm(formula=~team + direction),
                data = b_80[b_80$speed == 80 & b_80$Species == "B", ], meta.data = list(left = 50, width = 300))

region.table <- data.frame("Region.Label" = 1, "Area" = 117*0.3)
sample.table <- data.frame("Region.Label" = 1, "Sample.Label" = 1:16, "Effort" = 117)

dht(det_b_80, region.table = region.table, sample.table = sample.table[2, ], obs.table = b_80[b_80$Sample.Label == 2, ], options=list(convert.units=0.001))



