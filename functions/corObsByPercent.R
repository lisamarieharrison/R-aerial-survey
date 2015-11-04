corObsByPercent <- function(dat, percent_mat) {
  
  #corrects a sighting matrix by the % of sightings that observer 2 saw compared to observer 1
  #note = angle for inserted observations is the average angle for that species and trial
  #dat = matrix of sightings by observer 2
  #percent_mat = matrix of percentages of sightings compared to observer 1 with species as names
  #returns = matrix of corrected sightings
  
 
  for (t in unique(dat$Trial)) {
    for (s in c("B", "BOT", "S")) {
      cor_factor <- percent_mat[which(percent_mat[, 1] == s), 2]
      n_rows <- nrow(dat[dat$Trial == t & dat$Species == s, ])
      if (n_rows != 0) {
        n_rep <- abs(n_rows - round(n_rows / cor_factor))
        rep_row <- dat[dat$Trial == t & dat$Species == s, ][1, ]
        rep_row$Dist..from.transect <- mean(na.omit(dat[dat$Trial == t & dat$Species == s, ]$Dist..from.transect))
        for (i in 1:n_rep) {
          dat <- rbind(dat, rep_row)
        }
      }
    }
  }
  
  return(dat)
  
}
