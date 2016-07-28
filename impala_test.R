#RJ-MCMC for aerial survey data to run in parallel
#code modified from Cornelia Oedekoven
#author: Lisa-Marie Harrison
#date: 25/07/2016

runRjmcmc <- function (chain, init) {
  
  
  library(compiler)
  
  ########################################  set initial values ##############
  
  nind <- 73 #73 individuals
  nz <- 300
  
  # number of iterations
  nt <- 10000
  
  #initial parameters
  cur_theta <- 1
  cur_psi   <- 0.8
  
  #matrix to hold parameter estimates
  par <- matrix(NA, nrow = nt + 1, ncol = 4)
  colnames(par) <- c("theta", "phi", "N", "D")
  par[1, ] <- c(cur_theta, cur_phi, NA, NA)
  
  ############################### the priors ###########################################################
  
  l.prior.theta <- function (t) {
    
    theta.prior <- log(dunif(t, 0, 10))
    
    return (theta.prior)
  }
  
  l.prior.psi <- function (p) {
    
    psi.prior <- log(dunif(p, 0, 1))
    
    return (psi.prior)
    
  }
  
  ########################## the likelihood function ######################
  
  log.lik.fct <- function (p) {
    
    theta <- p[1] # theta
    psi   <- p[2] # psi
    
    #likelihood for the detection function
    fe <- log(sapply(dat$x[1:nind], f.hn.function, theta))
    
    #likelihood for the count model
    l.pois.bern <- log((psi ^ sum(dat$y))*(1 - psi)^(sum(dat$y == 0)))
    
    #return combined likelihood
    post <- sum(fe) + l.pois.bern
    return(post)
    
  }
  
  
  ########################## the detection function ######################
  
  # using a half-normal detection function for line transects
  f.hn.function <- function (dis, sigma) {
    f <- exp(-(dis/sigma)^2)
    return(f)
  }
  f.hn.function <- cmpfun(f.hn.function)
  
  #################################  the MCMC algorithm ######################################
  # nt is the number of iterations and is set above
  # row 1 is filled in with initial values for parameters and models
  
  for (i in 2:nt) {
    
    if (i %% 1000 == 0) {
      print(i)
    }
    
    ######### updating the density model parameters #############
    
    # theta
    new_theta <- cur_theta + rnorm(1, 0, 0.08) 
    num <- log.lik.fct(c(new_theta, cur_psi)) + l.prior.theta(new_theta) 
    den <- log.lik.fct(c(cur_theta, cur_psi)) + l.prior.theta(cur_theta) 
    A <- min(1, exp(num-den))
    if (runif(1) <= A) {
      cur_theta <- new_theta #if accepted, update current parameters
    } 
    
    ######### updating the count model parameters #############
    
    # psi
    new_psi <- cur_psi + max(rnorm(1, 0, 0.01), -cur_psi) # cannot become 0 or less
    num <- log.lik.fct(c(cur_theta, new_psi)) + l.prior.psi(new_psi)
    den <- log.lik.fct(c(cur_theta, cur_psi)) + l.prior.psi(cur_psi)
    A <- min(1, exp(num - den))
    if (runif(1) <= A) {
      cur_psi <- new_psi #if accepted, update current parameters
    } 
    
    ######### data augmentation #############
    
    z    <- rbinom(nind+nz, 1, cur_psi) # latent indicator variables from data augmentation
    x    <- runif(nind+nz, 0, 4) # distance is a random variable
    p <- f.hn.function(x, cur_theta)
    mu   <- z * p 
    dat$x <- x
    for (j in 1:length(mu)) {
      dat$y[j]   <- rbinom(1, 1, mu[j]) # observation model
    }
    
    N <- sum(z)
    D <- N / 48 # 48 km*km = total area of transects
    
    
    #save estimates
    par[i, ] <- c(cur_theta, cur_psi, N, D)
    
  } # end of iteration 
  
}

par <- na.omit(par)

summary_tab <- cbind(colMeans(par), apply(par, 2, sd))
rownames(summary_tab) <- c("theta", "psi", "N", "D")
colnames(summary_tab) <- c("mean", "sd")
summary_tab

par(mfrow = c(2, 2))

for (i in 1:4) {
  plot(par[, i])
}
