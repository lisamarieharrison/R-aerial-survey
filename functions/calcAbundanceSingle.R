calcAbundanceSingle <- function(area, dataframe, model) {
  
  #calculates density and abundance using distance sampling data
  #area = survey area (km2). Wollongong - Newcastle = 265km
  #dataframe = total_observations data frame from distanceData function
  #model = ddf model
  
  obs.table <- cbind(rep(1, nrow(dataframe)), dataframe$Trial, dataframe)
  colnames(obs.table)[1:2] <- c("Region.Label", "Sample.Label")
  
  region.table <- data.frame(matrix(c(1, area), ncol = 2, byrow = T))
  colnames(region.table) <- c("Region.Label", "Area")
  
  sample.table <- data.frame(cbind(rep(1, 54), c(1:54), rep(265, 54)))
  colnames(sample.table) <- c("Region.Label", "Sample.Label", "Effort")
  
  d <- dht(model, region.table = region.table, sample.table = sample.table, obs.table = obs.table, 
           options = list(convert.units = 0.001))
  
  return(d)
  
}