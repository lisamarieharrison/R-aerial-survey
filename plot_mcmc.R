#plotting mcmc results
#author: Lisa-Marie Harrison
#date: 25/07/2016

setwd("C:/Users/43439535/Documents/Lisa/phd/aerial survey/data")

model1 <- read.csv("models.csv", header = T)
model2 <- read.csv("models2.csv", header = T)

det_param1 <- read.csv("det.param.csv", header = T)
det_param2 <- read.csv("det.param2.csv", header = T)

count_param1 <- read.csv("count.param.csv", header = T)
count_param2 <- read.csv("count.param2.csv", header = T)


#detection model scale (depends on other covariates too)
plot(det_param1[, 1], type = "l", main = "Gamma detection function - scale", xlab = "iteration", ylab = "scale", ylim = c(100, 300))
points(det_param2[, 1], type = "l", col = "red")


#detection model shape 
plot(det_param1[, 2], type = "l", main = "Gamma detection function - shape", xlab = "iteration", ylab = "shape", ylim = c(5, 20))
points(det_param2[, 2], type = "l", col = "red")

#count model intercept
plot(count_param1[, 1], type= "l", main = "Poisson count model - intercept", xlab = "iteration", ylab = "intercept")
points(count_param2[, 1], type = "l", col = "red")

#detection model choice
det_models <- cbind(c(model1[1:100000, 1], model2[1:100000, 1]), rep(1:2, each = 100000))
table(det_models[, 1], det_models[, 2])

#count model choice
count_models <- cbind(c(model1[1:100000, 2], model2[1:100000, 2]), rep(1:2, each = 100000))
table(count_models[, 1], count_models[, 2])







