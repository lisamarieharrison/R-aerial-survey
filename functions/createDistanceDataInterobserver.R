createDistanceDataInterobserver <- function(species, overlap, obs1_missed, obs2_missed, truncate=NULL) {
  
  #sets up data frame for distance sampling analysis
  #species code:  1 = baitfish, 2 = bottlenose dolphin, 12 = sharks
  #truncate: optional truncation distance in m
  
  dataframe <- cbind(1:nrow(overlap), rep(1, nrow(overlap)), rep(1, nrow(overlap)), overlap$Dist..from.transect, overlap$Species, overlap$Trial, overlap$Number)
  dataframe <- rbind(dataframe, dataframe)
  dataframe[(nrow(overlap) + 1):nrow(dataframe), 2] <- 2
  
  dataframe <- rbind(dataframe, cbind((nrow(dataframe) + 1):(nrow(dataframe)  + nrow(obs1_missed)), rep(2, nrow(obs1_missed)), rep(1, nrow(obs1_missed)), obs1_missed$Dist..from.transect, obs1_missed$Species, obs1_missed$Trial, obs1_missed$Number))
  dataframe <- rbind(dataframe, cbind((nrow(dataframe) + 1):(nrow(dataframe)  + nrow(obs1_missed)), rep(1, nrow(obs1_missed)), rep(0, nrow(obs1_missed)), obs1_missed$Dist..from.transect, obs1_missed$Species, obs1_missed$Trial, obs1_missed$Number))
  dataframe <- rbind(dataframe, cbind((nrow(dataframe) + 1):(nrow(dataframe)  + nrow(obs2_missed)), rep(1, nrow(obs2_missed)), rep(1, nrow(obs2_missed)), obs2_missed$Dist..from.transect, obs2_missed$Species, obs2_missed$Trial, obs2_missed$Number))
  dataframe <- rbind(dataframe, cbind((nrow(dataframe) + 1):(nrow(dataframe)  + nrow(obs2_missed)), rep(2, nrow(obs2_missed)), rep(0, nrow(obs2_missed)), obs2_missed$Dist..from.transect, obs2_missed$Species, obs2_missed$Trial, obs2_missed$Number))
  
  dataframe[, 1] <- c(rep(1:nrow(overlap), 2), rep((nrow(overlap) + 1):((nrow(overlap) + 1)+nrow(obs1_missed) - 1), 2), rep((nrow(overlap) + 1 +nrow(obs1_missed)):(nrow(overlap) + nrow(obs1_missed) + nrow(obs2_missed)), 2))
  
  dataframe <- data.frame(dataframe)
  colnames(dataframe) <- c("object", "observer", "detected", "distance", "species", "Trial", "size")
  
  dataframe <- dataframe[dataframe$species == species, ]
  
  if (!is.null(truncate)) {
    dataframe <- dataframe[dataframe$distance <= truncate, ]
  }
  
  if (species != 2) {
    dataframe <- dataframe[, 1:6] #remove size unless species is Bottlenose dolphin
  }
  
  dataframe <- dataframe[dataframe$distance !=0 & !is.na(dataframe$distance), ]
  
  return(dataframe)
  
}