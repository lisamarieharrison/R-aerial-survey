correctLTRT <- function(dat) {
  
  #remove single "LT" or "RT" or times will be incorrect
  #dat = matrix of full data set
  
  dat <- dat[!dat$Time == "", ]
  dat[, c("Beaufort.Sea.State", "Water.clarity")] <- apply(dat[, c("Beaufort.Sea.State", "Water.clarity")], 2, as.integer)
  nrow_dat <- nrow(dat) - 1
  i = 2
  while(i < nrow(dat)) {
    if ((dat$Type[i] == "LT" & dat$Type[i - 1] == "LT") | (dat$Type[i] == "RT" & dat$Type[i - 1] == "RT")) {
      dat <- dat[-i, ]
      i <- i - 1
    }
    i <- i + 1
  }
  
  return(dat)
  
}
