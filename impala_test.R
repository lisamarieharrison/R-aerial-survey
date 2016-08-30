#RJ-MCMC test on impala data from Royale
#author: Lisa-Marie Harrison
#date: 25/07/2016

dat <- list(y=
              c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
            x=c(0.72, 0.26, 0.58, 0.92, 1.64, 0.85, 1.64, 1.57, 0.22, 0.72, 
                0.87, 0.51, 0, 0.73, 0, 1.29, 1.64, 0.72, 0.3, 0.71, 1.51, 0.69, 
                0.9, 0.65, 1.66, 0.38, 3.78, 0.78, 0.42, 0, 4, 1.75, 0.3, 0.35, 
                0.86, 0.32, 2, 2.72, 0.26, 0.77, 0.41, 2, 0.86, 0, 0.94, 0.55, 
                0.1, 0.85, 0, 0.78, 0, 0.96, 0, 0.64, 1.88, 0, 1.61, 1.5, 0.64, 
                1.93, 1.06, 1.15, 1.43, 1.29, 2.46, 1.23, 1.23, 1.53, 1.43, 0.34, 
                0.96, 2.6, 0.09, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA),
            nind=73,
            nz=300)

init <- list(theta = 1, psi = 0.8, z = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))



library(compiler)

########################################  set initial values ##############

nind <- 73 #73 individuals
nz <- 300

# number of iterations
nt <- 20000

#initial parameters
cur_theta <- init$theta
cur_psi   <- init$psi
z <- init$z
dat$x[(nind+1):length(dat$x)] <- runif(nz, 0, 4)

#matrix to hold parameter estimates
par <- matrix(NA, nrow = nt + 1, ncol = 4)
colnames(par) <- c("theta", "phi", "N", "D")
par[1, ] <- c(cur_theta, cur_psi, NA, NA)

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
  
  z <- rbinom(nind+nz, 1, psi) # latent indicator variables from data augmentation
  x <- runif(nind+nz, 0, 4)
  p     <- f.hn.function(x, theta)
  mu   <- z * p
  y <- NULL
  for (i in 1:length(mu)) {
    y[i] <- rbinom(1, 1, p[i])
  }
  
  #detection function likelihood 
  #eqn 5.18 (pg 66) in Distance Sampling: Methods and Applications
  u <- integrate(f.hn.function, 0, 4, theta)$value
  l.det <- prod(f.hn.function(dat$x[1:nind], theta)/u)

  #count likelihood 
  #l.count <- log(dbinom(sum(dat$y), sum(z), mean(mu[mu != 0])))
  
  loglik <- sum(l.det)
  if (is.nan(loglik) | is.infinite(loglik)) {
    loglik <- -100000
  }
  
  return(loglik)
  
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
  new_theta <- cur_theta + rnorm(1, 0, 0.05) 
  num <- log.lik.fct(c(new_theta, cur_psi)) + l.prior.theta(new_theta) 
  den <- log.lik.fct(c(cur_theta, cur_psi)) + l.prior.theta(cur_theta) 
  A <- min(1, exp(num-den))
  if (runif(1) <= A) {
    cur_theta <- new_theta #if accepted, update current parameters
  } 
  
  ######### updating the count model parameters #############
  
  # psi
  new_psi <- cur_psi + rnorm(1, 0, 0.01) # must be 0 - 1
  while (new_psi <=0 | new_psi >= 1) {
    new_psi <- cur_psi + rnorm(1, 0, 0.01) 
  }
  num <- log.lik.fct(c(cur_theta, new_psi)) + l.prior.psi(new_psi)
  den <- log.lik.fct(c(cur_theta, cur_psi)) + l.prior.psi(cur_psi)
  A <- min(1, exp(num - den))
  if (runif(1) <= A & num != -1e5) {
    cur_psi <- new_psi #if accepted, update current parameters
  } 
  
  ######### data augmentation #############
  
  z <- rbinom(nind+nz, 1, cur_psi) # latent indicator variables from data augmentation
  
  N <- sum(z)
  D <- N / 48 # 48 km*km = total area of transects
  
  
  #save estimates
  par[i, ] <- c(cur_theta, cur_psi, N, D)
  
} # end of iteration 



par <- na.omit(par)

par(mfrow = c(2, 2))
for (i in 1:4) {
  plot(par[, i])
}

#par <- par[1000:nrow(par), ] #remove 1000 sample burn-in

summary_tab <- cbind(colMeans(par), apply(par, 2, sd))
rownames(summary_tab) <- c("theta", "psi", "N", "D")
colnames(summary_tab) <- c("mean", "sd")
summary_tab
