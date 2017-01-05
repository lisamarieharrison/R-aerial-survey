#plots the three detection functions for theory chapter
#author: Lisa-Marie Harrison
#date: 05/01/2017

library(mrds)
library(fdrtool)

f.haz.function <- function(dis, sigma, shape) {
  f <- 1-exp(-(dis/sigma)^(-shape))
  return(f)
}


f.hn.function <- function(dis, sigma) {
  f <- exp(-dis^2/(2*sigma^2))
  f
}

f.gamma.function <- function (dis, key.scale, key.shape) {
  
  fr <- (1/gamma(key.shape)) * (((key.shape - 1)/exp(1))^(key.shape - 1))
  v1 <- dis/(key.scale * fr)
  return(v1^(key.shape-1)*exp(-v1)/(gamma(key.shape)*fr))
  
}

par(mfrow = c(1, 3), mar = c(2, 2, 2, 2), oma = c(4, 4, 1, 1))
plot(f.hn.function(1:1000, sigma = 1000), type = "l", main = "Half-normal", lwd = 2, xlab = "", ylab = "", cex.main = 2)
plot(f.haz.function(1:1000, sigma = 500, shape = 5), type = "l", main = "Hazard-rate", lwd = 2, xlab = "", ylab = "", cex.main = 2)
plot(f.gamma.function(1:1000, key.scale = 600, key.shape = 1.2), type = "l", main = "Gamma", lwd = 2, xlab = "", ylab = "", cex.main = 2)

mtext("Detection Probability", side = 2, outer = TRUE, line = 2, cex = 2)
mtext("Distance", side = 1, outer = TRUE, line = 2, cex = 2)


