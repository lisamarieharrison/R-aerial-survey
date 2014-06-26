setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/data")

data <- read.table(file = "sightings by envt.txt", header = F, sep = "\t", fill = T)

names(data) <- c("dist", "alt", "species", "number", "size", "swim_direc", 
                 "fish_yn", "dist_shore", "secondary", "sea_state", "wind_dir", 
                 "wind_speed", "cloud", "clarity", "other", "glare", "sst")

attach(data)

boxplot(sst ~ species, xaxt = "n", xlab = "Species", ylab = "SST (degrees)")
axis(side=1, at = 1:4, labels = c("Fish", "Bottlenose dolphin", "Humpback whale", 
                                  "White shark"))
title("Sea surface temperature by species")

barplot(table(species, sst), beside = T, xlab = "Temperature (degrees celcius)", ylab = "number of 
        sightings", col = c("blue", "red", "grey", "yellow"))
title("Number of sightings of each species by temperature")
legend("topleft", c("Fish,", "Bottlenose dolphin", "Humpback whale", "White shark"), bty = "n", 
       col = c("blue", "red", "grey", "yellow"), lwd = 2.5)