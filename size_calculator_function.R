#calculates the size of an animal using distance to animal and focal length

SizeCalc <- function(altitude, angle, focal_length, pixels, ccd_x, ccd_y){
    
  #calculate the distance to the animal from the aircraft
  air_distance    <- altitude/sin(angle*(pi/180))

  theta <- atan((ccd_x/2)/focal_length)
  d <- air_distance*tan(theta)
  pixel_size <- (2*d)/4948 #x direction number of pixels
  animal_size <- pixel_size * pixels
  
  return(animal_size)
  
}
