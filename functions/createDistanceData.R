createDistanceData <- function(species, lisa_obs, direction, truncate=NULL) {
  
  #creates data frame in the form mrds requires for distance sampling
  #species = character with the species code. e.g.: "BOT" = bottlenose dolphin
  #lisa_obs = data.frame of all observations
  #direction = character with flight direction. "N" = north, "S" = south
  #truncate = optional numeric specifying the distance (m) from transect to truncate sightings
  #return = data.frame with columns object, observer, detected, distance, species, Trial, size
  
  lisa_obs <- lisa_obs[lisa_obs$Flight.Direction == direction & 
                         lisa_obs$Type == "S" & 
                         lisa_obs$Dist..from.transect !=0 & 
                         !is.na(lisa_obs$Dist..from.transect) & 
                         lisa_obs$Species == species, ]
  
  if (!is.null(truncate)) {
    lisa_obs <- lisa_obs[lisa_obs$Dist..from.transect <= truncate, ]
  }
  
  total_observations <- cbind(1:nrow(lisa_obs), rep(1, nrow(lisa_obs)), rep(1, nrow(lisa_obs)), lisa_obs$Dist..from.transect, lisa_obs$Species, lisa_obs$Trial, lisa_obs$Number, lisa_obs$Season, lisa_obs$Year, lisa_obs$Beaufort.Sea.State, lisa_obs$Water.clarity, lisa_obs$Cloud.cover)
  total_observations <- data.frame(total_observations)
  colnames(total_observations) <- c("object", "observer", "detected", "distance", "species", "Trial", "size", "Season", "Year", "Sea.state", "Water.clarity", "Cloud.cover")
  
  total_observations <- apply(total_observations, 2, as.character)
  total_observations <- data.frame(apply(total_observations, 2, as.numeric))
  
  if (species == "BOT") {
    total_observations <- total_observations[!is.na(total_observations$size), ]
  } else {
    total_observations <- total_observations[, c(1:6, 7:ncol(total_observations))]
  }
  
  return(total_observations)
  
}
