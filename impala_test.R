#RJ-MCMC for aerial survey data to run in parallel
#code modified from Cornelia Oedekoven
#author: Lisa-Marie Harrison
#date: 25/07/2016

library(compiler)

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



########################################  set initial values ###########################

nind <- 73 #73 individuals
nz <- 300

# number of iterations
nt <- 20000

#initial parameters
theta <- init$theta
psi   <- init$psi
z <- init$z
dat$x[(nind+1):length(dat$x)] <- runif(nz, 0, 4)
x <- dat$x
y <- dat$y

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
  
  # sample theta
  dists <- NULL
  for (j in 1: length(seq(0, 10, by = 0.01))) {
    dists[j] <- mean(f.hn.function(x, seq(0, 10, by = 0.01)[j]))
  }
  theta <- seq(0, 10, by = 0.01)[which.max(dbinom(sum(dat$y), sum(z), dists*psi))] 

  # sample psi
  psi <- seq(0, 1, by = 0.01)[which.max(dbinom(sum(dat$y), sum(z), mean(f.hn.function(x, theta))*seq(0, 1, by = 0.01)))]

  ######### data augmentation #############
  
  z <- rbinom(nind+nz, 1, psi) # latent indicator variables from data augmentation
  y <- NULL
  x <- dat$x
  x[(nind+1):length(x)] <- runif(nz, 0, 4)
  p    <- f.hn.function(x, theta)
  mu   <- z * p
  for (j in 1:length(mu)) {
    y[j] <- rbinom(1, 1, mu[j])
  }
  
  N <- sum(z)
  D <- N / 48 # 48 km*km = total area of transects
  
  
  #save estimates
  par[i, ] <- c(theta, psi, N, D)
  
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

