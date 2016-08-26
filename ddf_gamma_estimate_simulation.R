#simulates distance data to validate how ddf estimates gamma detection function parameters
#author: Lisa-Marie Harrison
#date: 26/08/2016

set.seed(123)
gdata<-data.frame(object=1:1000,distance=rgamma(1000,scale=1,shape=1))

mm <- ddf(dsmodel=~cds(key="gamma"), data=gdata, method="ds",
          meta.data=list(width=max(gdata$distance)), control = list(parscale=FALSE, standardize=FALSE))

#check fitted parameters
mm



f.gamma.function <- function (distance, scale, shape) {
  
  fr <- (1/gamma(shape)) * (((shape - 1)/exp(1))^(shape - 1))
  v1 <- distance/(scale * fr)
  return(v1^(shape-1)*exp(-v1)/(gamma(shape)*fr))
  
}


plot(f.gamma.function(seq(0, max(gdata$distance), by = 0.1), 1, 1), type = "l")
points(f.gamma.function(0:max(total_observations$distance), 318.6604, 2.609934+1), type = "l", col = "blue")


#parameters just need to be exponentiated and then shape have 1 added