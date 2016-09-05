chain = 1
scale0 = 318
shape0 = 3.6
species = "BOT"
truncate_left = 0
truncate_right = 1000
int0 = 0.1
i = 1
nz <- 100

f.gamma.function <- function (distance, scale, shape) {
  
  fr <- (1/gamma(shape)) * (((shape - 1)/exp(1))^(shape - 1))
  v1 <- distance/(scale * fr)
  return(v1^(shape-1)*exp(-v1)/(gamma(shape)*fr))
  
}