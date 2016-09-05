

lik <- function(x) {
  
  z <- rbinom(nz, 1, x) # latent indicator variables from data augmentation
  dist <- runif(nz, 0, 1000)
  p     <- f.gamma.function(dist, sig1, sha2)
  mu   <- z * p
  
  #count likelihood 
  count <- NULL
  for (j in 1:length(counts)) {
    count[j] <- dbinom(counts[j], sum(z), mean(mu[mu != 0]))
  }
  
  return(log(prod(count)))
  
}

liks <- lapply(seq(0, 1, by = 0.001), lik)

plot(seq(0, 1, by = 0.001), liks)