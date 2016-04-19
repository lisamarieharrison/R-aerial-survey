setwd(dir = "~/Lisa/phd/aerial survey/data/speed trials/completed datasheets")
dat <- read.csv("Completed_matched_data.csv", header = TRUE)
library(MASS)
library(piecewiseSEM)

s <- "BOT"

count_lisa <- table(dat$Transect[dat$Lisa != 0 & dat$Species == s])
count_vic  <- table(dat$Transect[dat$Vic != 0 & dat$Species == s])

#fill in zeros when no animals were seen
count_lisa[which(!(1:15 %in% names(count_lisa)))] <- 0
count_lisa[which(1:15 %in% names(count_lisa))] <- table(dat$Transect[dat$Lisa != 0 & dat$Species == s])
count_vic[which(!(1:15 %in% names(count_vic)))] <- 0
count_vic[which(1:15 %in% names(count_vic))] <- table(dat$Transect[dat$Vic != 0 & dat$Species == s])

#speed for each trial by observer
lisa_speed <- c(80, 100, 100, 80, 80, 100, 100, 80, 80, 100, 100, 80, 80, 100, 100)
vic_speed  <- c(100, 80, 80, 100, 100, 80, 80, 100, 100, 80, 80, 100, 100, 80, 80)

#make data frame for models
counts <- c(count_lisa, count_vic)
speeds <- c(lisa_speed, vic_speed)
trial  <- rep(c(1:15), 2)
day <- rep(c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4), 2)

speed_dat <- cbind(counts, speeds, trial, day)
speed_dat <- speed_dat[order(trial, speeds), ] #order by trial and then speed
row.names(speed_dat) <- seq(1:nrow(speed_dat))
speed_dat <- as.data.frame(speed_dat)
speed_dat$counts <- as.integer(speed_dat$counts)
  
#before and after groups for t test
before <- speed_dat[seq(1, nrow(speed_dat), by = 2), 1]
after  <- speed_dat[seq(2, nrow(speed_dat), by = 2), 1]


#paired t test vs mixed model
t.test(before, after, paired = TRUE)

speed.lm <- lme(counts ~ speeds, random=~1|trial, data = speed_dat)
summary(speed.lm)

#mixed model with correlation
speed.lm <- lme(counts ~ speeds, random=~1|day/trial, correlation = corAR1(form=~1|day/trial), data = speed_dat)
summary(speed.lm)

#glmm
speed.lm <- glmm(counts ~ speeds, random=~1|trial, family.glmm = poisson.glmm, m = 100, varcomps.names = "trial", data = speed_dat)
summary(speed.lm)

#glmm with correlation structure
speed.lm <- glmmPQL(counts ~ speeds, random=~1|day/trial, family = poisson, data = speed_dat, correlation = corAR1(form=~1|day/trial), control = list(msMaxIter = 100))
summary(speed.lm)

#gamma glmm with correlation structure
speed.lm <- glmmPQL(counts ~ speeds, random=~1|day/trial, family = Gamma, data = speed_dat, correlation = corAR1(form=~1|day/trial), control = list(msMaxIter = 100))
summary(speed.lm)

#summary plots for glmmPQL
preds <- predict(speed.lm, type = "response")
plot(speed_dat$counts, preds)
plot(speed.lm, main = "residuals vs fitted")


#------------------------------- ALL SPECIES AT ONCE ----------------------------------#


count_lisa <- table(dat$Transect[dat$Lisa != 0], dat$Species[dat$Lisa != 0])[, c(1, 2, 7)]
count_vic  <- table(dat$Transect[dat$Vic != 0], dat$Species[dat$Vic != 0])[, c(1, 2, 7)]

#speed for each trial by observer
lisa_speed <- c(80, 100, 100, 80, 80, 100, 100, 80, 80, 100, 100, 80, 80, 100, 100)
vic_speed  <- c(100, 80, 80, 100, 100, 80, 80, 100, 100, 80, 80, 100, 100, 80, 80)

#make data frame for models
counts <- c(c(count_lisa), c(count_vic))
species <- rep(rep(c("B", "BOT", "S"), each = nrow(count_lisa)), 2)
speeds <- c(rep(lisa_speed, 3), rep(vic_speed, 3))
trial  <- rep(c(1:15), 6)
day <- rep(c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4), 6)
observer <- rep(c("Lisa", "Vic"), each = length(counts)/2)

speed_dat <- cbind(counts, speeds, trial, day, species, observer)
speed_dat <- speed_dat[order(trial, speeds), ] #order by trial and then speed
row.names(speed_dat) <- seq(1:nrow(speed_dat))
speed_dat <- as.data.frame(speed_dat)
speed_dat$counts <- as.integer(as.character(speed_dat$counts))
speed_dat$speeds <- as.numeric(as.character(speed_dat$speeds))
speed_dat$trial <- as.numeric(as.character(speed_dat$trial))

#before and after groups for t test
before <- speed_dat[speed_dat$speeds == 80 & speed_dat$species == "BOT", 1]
after  <- speed_dat[speed_dat$speeds == 100 & speed_dat$species == "BOT", 1]


#paired t test vs mixed model
t.test(before, after, paired = TRUE)


#mixed model with correlation
speed.lm <- lme(counts ~ species*speeds + observer, random=~1|day/trial, correlation = corAR1(form=~1|day/trial), data = speed_dat)
summary(speed.lm)


#mixed model with correlation
speed.lm <- lme(counts ~ speeds + observer, random=~1|day/trial, data = speed_dat[speed_dat$species == "B", ])
summary(speed.lm)

#r-squared
sem.model.fits(speed.lm, aicc = TRUE)

#glmm
speed.lm <- glmm(counts ~ speeds + observer, random=~1|trial, family.glmm = poisson.glmm, m = 100, varcomps.names = "trial", data = speed_dat[speed_dat$species == "B", ])
summary(speed.lm)


speed.lm <- glmmPQL(counts ~ speeds + observer, random=~1|trial, family = poisson, data = speed_dat[speed_dat$species == "B", ], correlation = corAR1(form=~1|trial), control = list(msMaxIter = 100))
summary(speed.lm)



