test_calcAbundanceSingle <- function() {
  
  library(mrds)
  source("C:/Users/Lisa/Documents/phd/aerial survey/R/code/R-aerial-survey/mrds_modified_functions.R")
  
  strip_width <- 1000
  total_observations <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/total_observations.csv", header = TRUE)
  
  p_total <- ddf(method = 'ds',dsmodel =~ cds(key = "gamma", formula=~1), 
                 data = total_observations, meta.data = list(left = 50, width = strip_width))
  suppressMessages(d <- calcAbundanceSingle(265*strip_width/1000, total_observations, p_total))
  
  checkEqualsNumeric(0.186, round(d$individuals$average.p, 3)) #check average.p
  checkEqualsNumeric(54, as.numeric(d$individuals$summary$k)) #check k
  checkEqualsNumeric(96.848, round(d$individuals$N$Estimate[1], 3)) #check abundance estimate

}

test_sightingsAtLevel <- function() {
  
  envt_var_south <- as.matrix(read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/envt_var_south.csv", header = TRUE))
  dat_south <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/dat_south.csv", header = TRUE)
  
  n_sightings <- sightingsAtLevel(environmenal_hrs = envt_var_south, dat = dat_south, env_variable = "Water.clarity", species = "B", flight_direction = "S")  
  
  checkEquals(c(16.7, 14.4, 13.0), round(c(n_sightings), 1))
  
}

test_createDistanceData <- function() {
  
  cor_obs <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/cor_obs.csv", header = TRUE)
  answer  <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/total_observations.csv", header = TRUE)
  
  question <- createDistanceData(species = "B", cor_obs, truncate = 1000, direction = "S")
  
  checkEquals(answer, question)
  
}

test_corObsByPercent <- function() {
  
  dat <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/vic_obs.csv", header = T)
  source("C:/Users/Lisa/Documents/phd/aerial survey/R/code/R-aerial-survey/percentOfLisa.R")
  percent_of_lisa <- percentOfLisa()  
  answer <- corObsByPercent(dat, percent) 
  
  checkEqualsNumeric(592, nrow(answer))
  
}

test_correctLTRT <- function() {
  
  dat_south <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/dat_south_ltrt.csv", header = T)
  
  answer <- correctLTRT(dat_south)
  
  checkEqualsNumeric(2908, nrow(answer))
  
}

test_calcEnvtEffort <- function() {
  
  dat_south <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/dat_south.csv", header = TRUE)
  
  suppressMessages(envt.var.south <- calcEnvtEffort(dat_south))
  
  checkEqualsNumeric(4216, sum(envt.var.south[1:4]))
  
}

test_checkSameSighting <- function() {
  
  dat <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/vic_obs.csv", header = T)
  
  checkTrue(checkSameSighting(dat[2, ], dat[2, ]))
  checkTrue(!checkSameSighting(dat[2, ], dat[492, ]))
  
}

test_findOverlappingObservations <- function() {
  
  obs1 <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/lisa_obs_subset.csv", header = T)
  obs2 <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/vic_obs_subset.csv", header = T)
  
  obs <- findOverlappingObservations(set1 = obs1, set2 = obs2)
    
  checkEqualsNumeric(8, nrow(obs$overlap))
  checkEqualsNumeric(3, nrow(obs$obs1_missed))
  checkEqualsNumeric(16, nrow(obs$obs2_missed))
  
}

test_createDistanceDataInterobserver <- function() {
  
  obs1 <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/lisa_obs_missed_subset.csv", header = T)
  obs2 <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/vic_obs_missed_subset.csv", header = T)
  overlap <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/overlap_sub.csv", header = T)
  
  total_observations <- createDistanceDataInterobserver(1, overlap, obs1, obs2, truncate = 1000)
  
  checkEqualsNumeric(28, nrow(total_observations))
  
}

test_calcAbundanceInterobserver <- function() {
  
  library(mrds)
  
  total_observations <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/total_observations_interobserver.csv", header = T)
  
  p_total <- ddf(method="io", mrmodel =~ glm(~distance), dsmodel =~ cds(key = "gamma", formula=~1),
                 data = total_observations, meta.data = list(left = 50, width = 1000, point = FALSE))
  
  suppressMessages(d <- calcAbundanceInterobserver(265, dataframe = total_observations, p_total, type = p_total$method))
  
  checkEqualsNumeric(158, round(d$individuals$N$Estimate))
  checkEqualsNumeric(491, d$individuals$summary$n)
  checkEqualsNumeric(19, d$individuals$summary$k)
  
}

test_poolAllSharks <- function() {
  
  dat <- read.csv("C:/Users/Lisa/Documents/phd/aerial survey/R/data/test_data/all_sharks.csv", header = T)
  
  checkTrue("W" %in% levels(dat$Species))
  
  dat <- poolAllSharks(dat)
  
  checkTrue(!"W" %in% levels(dat$Species))
  
}






