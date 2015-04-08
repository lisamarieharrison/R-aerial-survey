#inter-observer aerial survey preliminary analysis
setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/R/data")
dat <- read.csv("interobserver_20150407.csv", header = T)

#total sightings seen by each observer
table(dat$Observer, dat$Species)

#number of overlapping sightings
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



findOverlap <- function(set1, set2, i, j, overlap, lisa_missed, vic_missed) {
  

  set2 <- set2[set2$Date == set1$Date[i], ]
  
  if (j > nrow(set2)) {
    print(paste(i, "No overlapping observation found. Try next observation"))
    lisa_missed <- rbind(lisa_missed, c(as.character(set1$Date[i]), as.character(set1$Species[i])))
    return(list(overlap = overlap, vic_missed = vic_missed, lisa_missed = lisa_missed))
  }
  
  if (checkSame(set1[i, ], set2[j, ])) {
    rbind(overlap, set1[i, ])
    print("Success: found a same observation!")
    vic_missed <- vic_missed[-which(vic_missed$Date == set1$Date[i] & substr(vic_missed$Time, 1, 5) == substr(set1$Time[i], 1, 5) & vic_missed$Species == set1$Species[i]), ]
    return(list(overlap = overlap, vic_missed = vic_missed, lisa_missed = lisa_missed))
  } else {
    findOverlap(set1, set2, i = i, j = j + 1)
  }
    
}

overlap <- dat[0, ]
vic_missed <- lisa
lisa_missed <- matrix(0, ncol = 2, nrow = 1)

for (i in 1:nrow(vic)) {
  result <- findOverlap(set1 = vic, set2 = vic_missed, i = i, j = 1, 
                        overlap = overlap, lisa_missed = lisa_missed, 
                        vic_missed = vic_missed)
  overlap <- result$overlap
  lisa_missed <- result$lisa_missed
  vic_missed <- result$vic_missed
}





