#calculates the size of an animal
setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/R")

altitude      <- 33 
angle         <- 55
focal_length 	<- 18
pixels        <- 377.2
ccd_x         <- 23.6
ccd_y         <- 15.6

source("size_calculator_function.R")


SizeCalc(altitude, angle, focal_length, pixels, ccd_x, ccd_y)


