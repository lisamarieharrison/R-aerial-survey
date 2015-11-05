findOverlappingObservations <- function(set1, set2) {
  
  #find overlapping observations and unique observations from 2 observers
  #set1 = full data frame of observations from observer 1
  #set2 = full data frame of observations from observer 2
  
  overlap <- dat[0, ]
  obs2_missed <- dat[0, ]
  obs1_missed <- dat[0, ]
  
  #remove observations with no time
  set1 <- set1[!set1$Time == "", ]
  set2 <- set2[!set2$Time == "", ]
  
  for (i in unique(set1$Date)) {
    
    obs1_temp <- set1[set1$Date == i, ]
    obs2_temp <- set2[set2$Date == i, ]

    for (r in 1:nrow(obs1_temp)) {
      
      r2 <- 1
      
      #check for match in obs2_temp
      while (r2 <= nrow(obs2_temp)) {
        
        if (checkSameSighting(obs1_temp[r, ], obs2_temp[r2, ])) { #if same
          
          overlap <- rbind(overlap, obs1_temp[r, ]) #add row to overlap
          obs2_temp <- obs2_temp[-r2, ] #remove row from obs2_temp
          r2 <- nrow(obs2_temp) + 1 #force break
          
        } else { #if not same
          if (r2 == nrow(obs2_temp)) {
            obs2_missed <- rbind(obs2_missed, obs1_temp[r, ]) #add to obs2_missed
          }
          r2 <- r2 + 1
        }
                  
      }
      
    }
    
    #add remaining vic sightings to obs1_missed
    obs1_missed <- rbind(obs1_missed, obs2_temp)
    
  }
  
  return(list(overlap = overlap, obs1_missed = obs1_missed, obs2_missed = obs2_missed))

}
