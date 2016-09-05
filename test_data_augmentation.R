

lik <- function(x) {
  
  z <- rbinom(nz, 1, x) # latent indicator variables from data augmentation
  dist <- runif(nz, 0, 1000)
  p     <- f.gamma.function(dist, scale0, shape0)
  mu   <- z * p
  
  #count likelihood 
  count<- dbinom(round(mean(counts)), sum(z), mean(mu[mu != 0]))

  
  return(log(count))
  
}

liks <- lapply(seq(0, 1, by = 0.001), lik)

plot(seq(0, 1, by = 0.001), liks)