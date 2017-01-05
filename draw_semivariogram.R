#Diagram of a semivariogram for theory chapter
#author: Lisa-Marie Harrison
#date: 05/01/2017


dist <- 1:100


par(mar = c(5, 5, 2, 2))

plot(dist, (2*dist + 1)/(dist + 5), type = "l", ylim = c(0, 2.5), lwd = 2, bty = "l", xlab = "Lag distance", ylab = "Semivariance", cex.lab = 2, xaxt = "n", yaxt = "n")
points(c(-2, 80), c(1.9, 1.9), type = "l", lty = 2)
points(c(80, 80), c(-2, 1.9), type = "l", lty = 2)
points(c(-2, 100), c(0.5, 0.5), type = "l", lty = 2)

text(93, 1.3, "Range", cex = 2)
text(20, 0.1, "Nugget", cex = 2)
text(20, 2.35, "Sill", cex = 2)

arrows(20, 0.2, 20, 0.4, length = 0.1, lwd = 2)
arrows(20, 2.2, 20, 1.95, length = 0.1, lwd = 2)
arrows(91, 1.2, 81, 1, length = 0.1, lwd = 2)