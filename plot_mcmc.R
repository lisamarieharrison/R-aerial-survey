#plotting mcmc results
#author: Lisa-Marie Harrison
#date: 25/07/2016

setwd("C:/Users/43439535/Documents/Lisa/phd/aerial survey/data")
library(mcmcse)

model1 <- read.csv("models1.csv", header = T)
model2 <- read.csv("models2.csv", header = T)

det_param1 <- read.csv("det.param1.csv", header = T)
det_param2 <- read.csv("det.param2.csv", header = T)

count_param1 <- read.csv("count.param1.csv", header = T)
count_param2 <- read.csv("count.param2.csv", header = T)

N1 <- read.csv("N1.csv", header = F)
N2 <- read.csv("N2.csv", header = F)

#detection model scale (depends on other covariates too)
det_param1 <- na.omit(det_param1)
det_param2 <- na.omit(det_param2)
plot(det_param1[, 1], type = "l", main = "Gamma detection function - scale", xlab = "iteration", ylab = "scale", ylim = c(min(det_param1[, 1], det_param2[, 1]), max(det_param1[, 1], det_param2[, 1])))
points(det_param2[, 1], type = "l", col = "red")


#detection model shape 
plot(det_param1[, 2], type = "l", main = "Gamma detection function - shape", xlab = "iteration", ylab = "shape", ylim = c(min(det_param1[, 2], det_param2[, 2]), max(det_param1[, 2], det_param2[, 2])))
points(det_param2[, 2], type = "l", col = "red")

#count model intercept
count_param1 <- na.omit(count_param1)
count_param2 <- na.omit(count_param2)
plot(count_param1[, 1], type= "l", main = "Count model - intercept", xlab = "iteration", ylab = "intercept", ylim = c(min(count_param1[, 1], count_param2[, 1]), max(count_param1[, 1], count_param2[, 1])))
points(count_param2[, 1], type = "l", col = "red")

#abundance estimate
plot(N1$V1, type = "l", ylim = c(0, 150))
points(N2$V1, type = "l", col = "red")


#detection model choice
det_models <- cbind(c(model1[1:100000, 1], model2[1:100000, 1]), rep(1:2, each = 100000))
table(det_models[, 1], det_models[, 2])

#count model choice
count_models <- cbind(c(model1[1:100000, 2], model2[1:100000, 2]), rep(1:2, each = 100000))
table(count_models[, 1], count_models[, 2])

#choose burn in 
burn_in <- 1000




#---------------- CHECK PARAMETER CONVERGENCE ---------------------#

gelman.rubin <- function(param) {
  # mcmc information
  n <- nrow(param) # number of iterations
  m <- ncol(param) # number of chains
  
  # calculate the mean of the means
  theta.bar.bar <- mean(colMeans(param))
  
  # within chain variance
  W <- mean(apply(param, 2, var))
  
  # between chain variance
  B <- n / (m - 1) * sum((colMeans(param) - theta.bar.bar) ^ 2)
  
  # variance of stationary distribution
  theta.var.hat <- (1 - 1 / n) * W + 1 / n * B
  
  # Potential Scale Reduction Factor (PSRF)
  R.hat <- sqrt(theta.var.hat / W)
  
  return(R.hat)
}

#parameter considered converged if Gelman-Rubin statistic is < 1.1
gelman.rubin(cbind(det_param1[1:burn_in, 1], det_param2[1:burn_in, 1]))


#--------------------- CHECK DETECTION FUNCTION FIT -----------------------#

#plot detection function over histogram of distances
#distances scaled by expected count
#remove burn in and average over all other parameter values (doesn't take model choice into account)



if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "~/Lisa/phd/aerial survey/data")
  source_location <- "~/Lisa/phd/aerial survey/R/R-aerial-survey/functions/"
} else {
  setwd(dir = "~/phd/aerial survey/data")
  source_location <- "~/Documents/phd/aerial survey/R/R-aerial-survey/functions/"
}

source(paste0(source_location, "gamma_det_fun.R"))

dat <- read.csv("aerial_survey_summary_r.csv", header = T) #lisa's sighting data
dat$Year <- as.numeric(substr(as.character(dat$Date), nchar(as.character(dat$Date)) - 3, nchar(as.character(dat$Date))))
covey.d <- dat[dat$Species == "BOT" & !is.na(dat$Dist..from.transect) & dat$Secondary != "Y" & dat$Observer == "Lisa", c("Date", "Season", "Dist..from.transect", "Beaufort.Sea.State", "Cloud.cover", "Water.clarity", "Flight.Direction")]
covey.d$Date <- as.character(covey.d$Date)
covey.d$Year <- as.numeric(substr(covey.d$Date, nchar(covey.d$Date) - 3, nchar(covey.d$Date))) #extract year from date
names(covey.d) <- c("Date", "Season", "Distance", "Sea_state", "Cloud_cover", "Water_clarity", "Flight.Direction", "Year")
covey.d$Id <- 1:nrow(covey.d)
covey.d$Visit <- as.numeric(as.factor(covey.d$Date)) #visit is transect
covey.d <- covey.d[covey.d$Distance != 0, ] #remove 0 distances because they are errors
covey.d <- covey.d[covey.d$Distance <= 1000 & covey.d$Flight.Direction == "S", ] #truncate to 1km

det.param  <- as.matrix(na.omit(det_param1))
det.param2 <- as.matrix(na.omit(det_param2))

hist.obj <- hist(covey.d$Distance, plot = FALSE, breaks = 15)

#calculate scale averaged across all parameter levels
calc_scale1 <- det.param[burn_in:nrow(det.param), 1] * exp(sum(colMeans(det.param[burn_in:nrow(det.param), 3:18])))

calc_scale2 <- det.param2[burn_in:nrow(det.param2), 1] * exp(sum(colMeans(det.param2[burn_in:nrow(det.param2), 3:18])))

scale <- mcse(c(calc_scale1, calc_scale2), size = "sqroot", g = NULL,
              method = c("bm", "obm", "tukey", "bartlett"),
              warn = FALSE)

shape <- mcse(c(det.param[burn_in:nrow(det.param), 2], det.param2[burn_in:nrow(det.param2), 2]), size = "sqroot", g = NULL,
              method = c("bm", "obm", "tukey", "bartlett"),
              warn = FALSE)

scale_mean <- scale$est
shape_mean <- shape$est

nc <- length(hist.obj$mids)
pa <- integrate(f.gamma.function, 0, max(covey.d$Distance), scale=scale_mean, shape=shape_mean)$value/max(covey.d$Distance)
Nhat <- nrow(covey.d)/pa
breaks <- hist.obj$breaks
expected.counts <- (breaks[2:(nc+1)]-breaks[1:nc])*(Nhat/breaks[nc+1])


nc <- length(breaks)-1
pdot <- f.gamma.function(covey.d$Distance, scale=scale_mean, shape=shape_mean)
Nhat <- sum(1/pdot)


hist.obj$density <- hist.obj$counts/expected.counts
hist.obj$density[expected.counts==0] <- 0
hist.obj$equidist <- FALSE

plot(hist.obj, ylim = c(0, 1.3), xlim = c(0, 1000), col = "lightgrey", xlab = "Distance (m)")
points(f.gamma.function(0:1000, scale=scale_mean, shape=shape_mean), type = "l", col = "red", lwd = 2)


#----------------------- PARAMETER ESTIMATES --------------------------#

scale_se <- scale$se
shape_se <- shape$se

scale_cv <- scale_se/scale_mean
shape_cv <- shape_se/shape_mean

#summary table
sum_sats <- rbind(c(shape_mean, shape_se, shape_cv), c(scale_mean, scale_se, scale_cv))
colnames(sum_sats) <- c("Est", "SD", "CV")
rownames(sum_sats) <- c("Shape", "Scale")
sum_sats <- round(sum_sats, 4)

sum_sats


#average p
average.p <- integrate(f.gamma.function, 0, 1000, scale_mean, shape_mean)$value/1000
average.p

#calculate abundance
#currently doesn't support other coefficients

#data augmentation

Ns <- na.omit(c(N1[burn_in:nrow(N1), ], N2[burn_in:nrow(N2), ]))

n_hat <- mean(Ns)
n_sd <- sd(Ns)/sqrt(length(Ns))
n_cv <- n_sd/n_hat

paste0("N = ", round(n_hat, 3), " (CV = ", round(n_cv, 3), ")")


