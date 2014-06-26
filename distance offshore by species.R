setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/R")

data <- read.csv(file = "season aerial survey.csv", header = T)
attach(data)

#bin data into 25m intervals
bin <- cut(data[, 13], breaks = seq(0, 550, 25))
df <- table(data[, 15], bin)

plot(as.numeric(df[1, ]), xaxt = "n", xlab = "Distance offshore (m)",
     ylab = "Number of Sightings", type = "l", lwd = 2)
points(as.numeric(df[2, ]), type = "l", col = "red", lwd = 2)
points(as.numeric(df[3, ]), type = "l", col = "blue", lwd = 2)
axis(side = 1, labels = levels(bin), at = seq(1, 22))
title("Sightings by distance offshore for southbound transects")

legend("topright", c("Fish", "Bottlenose dolphin", "Humpback whale"), 
       col = c("black", "red", "blue"), lwd = 2.5, bty = "n")


boxplot(Dist..from.transect ~ Species, xlab = "Species", ylab = "Distance offshore (m)", xaxt = 'n')
axis(side = 1, at = 1:3, c("Fish", "Bottlenose dolphin", "Humpback whale"))
title("Distance offshore by species for Southbound transects")




