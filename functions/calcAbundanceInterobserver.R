calcAbundanceInterobserver <- function(area, dataframe, model, type) {
  
  #calculates density and abundance using distance sampling data
  #area = survey area (km2). Wollongong - Newcastle = 265km
  #dataframe = total_observations data frame from distanceData function
  #model = ddf model
  #type = "character containing which obs to include, "trial" = only observer 1, "io" = both
  
  
  obs.table <- cbind(rep(1, nrow(dataframe)), dataframe$Trial, dataframe)
  colnames(obs.table)[1:2] <- c("Region.Label", "Sample.Label")
  
  if (type == "trial") {
    obs.table <- obs.table[obs.table$observer == 1 & obs.table$detected == 1, ]
  } else {
    obs.table <- obs.table[duplicated(obs.table$object), ] #remove duplicate detections
  }
  
  region.table <- data.frame(matrix(c(1, area), ncol = 2, byrow = T))
  colnames(region.table) <- c("Region.Label", "Area")
  
  sample.table <- data.frame(cbind(rep(1, length(unique(dataframe$Trial))), c(1:length(unique(dataframe$Trial))), rep(265, length(unique(dataframe$Trial)))))
  colnames(sample.table) <- c("Region.Label", "Sample.Label", "Effort")
  
  d <- dht(model, region.table = region.table, sample.table = sample.table, obs.table = obs.table, 
           options = list(convert.units = 0.001))
  
  return(d)
  
}