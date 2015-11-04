sightingsAtLevel <- function(environmenal_hrs, dat, env_variable, species, flight_direction) {
  
  #calculates the number of sightings per hour at each environmental condition level
  #environmental_hrs = vector of hours at each level. e.g.: envt.var.south
  #dat = full data set of observations. e.g.: lisa_obs
  #env_variable = character containing the name of the environmental variable
  #species = character containing the species code. e.g.: "B" = baitfish
  #flight_direction = character containing the flight direction. "S" = south, "N" = north
  
  hrs <- environmenal_hrs[grep(env_variable, colnames(environmenal_hrs))]/60
  t <- table(factor(dat[dat$Species == species & dat$Flight.Direction == flight_direction, which(colnames(dat) == env_variable)],
                    levels = unique(dat[dat$Flight.Direction == flight_direction, 
                                        which(colnames(dat) == env_variable)])))
  t <- t[order(names(t))]                                                                                                                                 
  sightings_hr <- matrix(t/hrs, ncol = length(t))
  colnames(sightings_hr) <- names(t)
  
  #species hash table
  species_full <- cbind(c("B", "BOT", "S"), c("Baitfish", "Dolphins", "Sharks"))
  Species <- species_full[which(species_full[, 1] == species), 2]
  
  return(sightings_hr)
  
}
