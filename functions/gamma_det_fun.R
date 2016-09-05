#gamma detection function for distance sampling analysis
#author: Lisa-Marie Harrison
#date: 5/9/2016
#modified from mrds package

f.gamma.function <- function (distance, scale, shape) {
  
  fr <- (1/gamma(shape)) * (((shape - 1)/exp(1))^(shape - 1))
  v1 <- distance/(scale * fr)
  return(v1^(shape-1)*exp(-v1)/(gamma(shape)*fr))
  
}