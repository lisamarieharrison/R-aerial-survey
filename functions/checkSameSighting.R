checkSameSighting <- function(r1, r2) {
  #checks if 2 observations are the same using time and date
  #r1, r2 = row of single observation from the full observation matrix
  #returns: boolean
  
  if (r1$Date == r2$Date & substr(r1$Time, 1, 5) == substr(r2$Time, 1, 5) & as.character(r1$Species) == as.character(r2$Species)) {
    return(TRUE)
  } 
  
  #add extra leniency if sightings are both dolphins and are within 1 min of each other because they can be spread out
  if (r1$Species == "BOT" & r2$Species == "BOT" & abs(as.numeric(chron(times. = as.character(r1$Time), format = "h:m:s") - chron(times. = as.character(r2$Time), format = "h:m:s")))*24*60 < 1) {
    return(TRUE)
  }  else {
    return(FALSE)
  }
}
  