#plotting mcmc results
#author: Lisa-Marie Harrison
#date: 25/07/2016

setwd("C:/Users/43439535/Documents/Lisa/phd/aerial survey/data")

model1 <- read.csv("models1.csv", header = T)
model2 <- read.csv("models2.csv", header = T)

det_param1 <- read.csv("det.param1.csv", header = T)
det_param2 <- read.csv("det.param2.csv", header = T)

count_param1 <- read.csv("count.param1.csv", header = T)
count_param2 <- read.csv("count.param2.csv", header = T)


#detection model scale (depends on other covariates too)
det_param1 <- na.omit(det_param1)
det_param2 <- na.omit(det_param2)
plot(det_param1[, 1], type = "l", main = "Gamma detection function - scale", xlab = "iteration", ylab = "scale", ylim = c(min(det_param1[, 1], det_param2[, 1]), max(det_param1[, 1], det_param2[, 1])))
points(det_param2[, 1], type = "l", col = "red")


#detection model shape 
plot(det_param1[, 2], type = "l", main = "Gamma detection function - shape", xlab = "iteration", ylab = "shape", ylim = c(min(det_param1[, 2], det_param2[, 2]), max(det_param1[, 2], det_param2[, 2])))
points(det_param2[, 2], type = "l", col = "red")

#count model intercept
plot(count_param1[, 1], type= "l", main = "Poisson count model - intercept", xlab = "iteration", ylab = "intercept", ylim = c(0, 5))
points(count_param2[, 1], type = "l", col = "red")


#detection model choice
det_models <- cbind(c(model1[1:100000, 1], model2[1:100000, 1]), rep(1:2, each = 100000))
table(det_models[, 1], det_models[, 2])

#count model choice
count_models <- cbind(c(model1[1:100000, 2], model2[1:100000, 2]), rep(1:2, each = 100000))
table(count_models[, 1], count_models[, 2])

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
gelman.rubin(cbind(det_param1[1:10000, 1], det_param2[1:10000, 1]))


#--------------------- CHECK DETECTION FUNCTION FIT ------------------#

#plot detection function over histogram of distances
#distances scaled by expected count



if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/aerial survey/data")
  source_location <- "~/Lisa/phd/Mixed models/R code/R-functions-southern-ocean/"
} else {
  setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/data")
  source_location <- "~/phd/southern ocean/Mixed models/R code/R-functions-southern-ocean/"
}

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


f.gamma.function <- function (distance, scale, shape) {
  
  fr <- (1/gamma(shape)) * (((shape - 1)/exp(1))^(shape - 1))
  v1 <- distance/(scale * fr)
  return(v1^(shape-1)*exp(-v1)/(gamma(shape)*fr))
  
}

det.param <- as.matrix(na.omit(det_param1))

hist.obj <- hist(covey.d$Distance, plot = FALSE)

#calculate scale averaged across all parameter levels
calc_scale <- det.param[nrow(det.param) - 1, 1] * exp(mean(c(0, det.param[, 3:4])) +
                                                        mean(c(0, det.param[, 5:6])) +
                                                        mean(c(0, det.param[, 7:8])) +
                                                        mean(c(0, det.param[, 9:16])) +
                                                        mean(c(0, det.param[, 17:18])))

shape <- det.param[nrow(det.param) - 1, 2]

nc <- length(hist.obj$mids)
pa <- integrate(f.gamma.function, 0, max(covey.d$Distance), scale=calc_scale, shape=det.param[nrow(det.param) - 1, 2])$value/max(covey.d$Distance)
Nhat <- nrow(covey.d)/pa
breaks <- hist.obj$breaks
expected.counts <- (breaks[2:(nc+1)]-breaks[1:nc])*(Nhat/breaks[nc+1])


nc <- length(breaks)-1
pdot <- f.gamma.function(covey.d$Distance, scale=calc_scale, shape=shape)
Nhat <- sum(1/pdot)



hist.obj$density <- hist.obj$counts/expected.counts
hist.obj$density[expected.counts==0] <- 0
hist.obj$equidist <- FALSE

plot(hist.obj, ylim = c(0, 1))
points(f.gamma.function(0:max(covey.d$Distance), scale=calc_scale, shape=shape), type = "l", col = "red")

