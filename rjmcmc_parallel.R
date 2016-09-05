#RJ-MCMC for aerial survey data to run in parallel
#code modified from Cornelia Oedekoven
#author: Lisa-Marie Harrison
#date: 25/07/2016

runRjmcmc <- function (chain, scale0, shape0, int0, species, nz, truncate_left, truncate_right) {
  
  library(compiler)
  library(tcltk)
  
  if (Sys.info()[4] == "SCI-6246") {
    setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/aerial survey/data")
    source_location <- "~/Lisa/phd/Mixed models/R code/R-functions-southern-ocean/"
  } else {
    setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/data")
    source_location <- "~/phd/southern ocean/Mixed models/R code/R-functions-southern-ocean/"
  }
  
  dat <- read.csv("aerial_survey_summary_r.csv", header = T) #lisa's sighting data
  dat$Year <- as.numeric(substr(as.character(dat$Date), nchar(as.character(dat$Date)) - 3, nchar(as.character(dat$Date))))
  covey.d <- dat[dat$Species == species & !is.na(dat$Dist..from.transect) & dat$Secondary != "Y" & dat$Observer == "Lisa", c("Date", "Season", "Dist..from.transect", "Beaufort.Sea.State", "Cloud.cover", "Water.clarity", "Flight.Direction")]
  covey.d$Date <- as.character(covey.d$Date)
  covey.d$Year <- as.numeric(substr(covey.d$Date, nchar(covey.d$Date) - 3, nchar(covey.d$Date))) #extract year from date
  names(covey.d) <- c("Date", "Season", "Distance", "Sea_state", "Cloud_cover", "Water_clarity", "Flight.Direction", "Year")
  covey.d$Id <- 1:nrow(covey.d)
  covey.d$Visit <- as.numeric(as.factor(covey.d$Date)) #visit is transect
  covey.d <- covey.d[covey.d$Distance != 0, ] #remove 0 distances because they are errors
  covey.d <- covey.d[covey.d$Distance <= truncate_right & covey.d$Distance >= truncate_left & covey.d$Flight.Direction == "S" & covey.d$Sea_state < 4, ] #truncate to 1km
  
  length.d <- length(covey.d$Distance)
  
  
  cd <- data.matrix(covey.d) #required for apply functions
  
  visit_tab <- data.matrix(covey.d[match(1:max(covey.d$Visit), covey.d$Visit), ])
  
  for (day in which(is.na(visit_tab[, 10]))) {
    
    visit_tab[day, ] <- c(NA, unique(dat$Season[as.numeric(dat$Date) == day]), NA, NA, NA, NA, NA, unique(dat$Year[as.numeric(dat$Date) == day]),  NA, day)
    
  }
  counts <- rep(0, max(covey.d$Visit))
  counts[unique(covey.d$Visit)] <- table(covey.d$Visit)
  
  rm(dat)
  
  years          <- sort(unique(covey.d$Year))
  seasons        <- levels(covey.d$Season)

  j <- length(counts)
  
  ########################################  set initial values ##############
  
  #  count model parameters: intercept and random effect standard deviation for L_n(\bmath{\beta}|\bmath{\theta}) (eqn (6))
  std.ran0 <- 1
  
  # the random effect coefficients b_j
  b0 <- rnorm(j, 0, std.ran0)
  
  #########################################################################
  # setting up the matrices that will contain the paramter values;
  
  # number of iterations
  nt <- 20000 
  
  df_size <- 100
  
  # holds the values for detection function parameters for each iteration
  # for det model we have 23 parameters: scale, shape, year(3), season(3), sea_state(3), cloud(9), water_clarity(3)
  # for count model we have 7 parameters: intercept, year(3) and season(3)
  det.param <- matrix(NA, df_size + 1, 18)
  
  # for an intercept only model:
  det.param[2, ] <- c(scale0, shape0, rep(0, 16))
  
  # holds the model id number for detection function for each iteration, 
  det.model <- matrix(NA, df_size + 1, 1)            # refers to det.list
  
  # the matrix that will keep the parameter values for the density model
  count.param <- matrix(NA, df_size + 1, 8 + j)
  
  #matrix to hold abundance estimates
  N <- matrix(NA, df_size + 1, 1)
  
  # filling in the initial values
  count.param[2, ] <- c(int0, rep(0, 6), std.ran0, b0)
  
  # holds the model id number for density model for each iteration
  count.model <- matrix(NA, df_size + 1, 1)          # refers to count.list
  
  ############# proposal distributions
  # proposal distributions for detection function parameters:
  
  msyt.prop.mean <- c(320, 4, rep(0.5, 16))
  msyt.prop.sd <- c(10, rep(0.5, 17))
  
  
  # model identifier for detection function, model number in rows, parameters (y/n) in columns
  det.list <- matrix(NA, 9, 18)
  colnames(det.list) <- c("sig","sha", "year2013", "year2014", "seasonSummer", "seasonSpring", paste0("s", 1:2),  paste0("cc", 0:7), paste0("wc", 1:2))
  det.list[1, ] <-c(rep(1, 4), rep(0, 14))                # mcds with year 
  det.list[2, ] <-c(1, 1, 0, 0, 1, 1, rep(0, 12))         # mcds with season
  det.list[3, ] <-c(1, 1, rep(0, 4), 1, 1, rep(0, 10))    # mcds with sea state
  det.list[4, ] <-c(1, 1, rep(0, 14), 1, 1)               # mcds with water clarity
  det.list[5, ] <-c(1, 1, rep(0, 6), rep(1, 8), 0, 0)     # mcds with cloud cover
  det.list[6, ] <-c(rep(1, 6), rep(0, 12))                # mcds with year and season
  det.list[7, ] <-c(rep(1, 8), rep(0, 10))                # mcds with year, season and sea state
  det.list[8, ] <-c(1, 1, rep(0, 16))                     # global hazard rate model
  det.list[9, ] <-c(1, 1, rep(1, 16))                     # mcds with all variables
  
  
  # model identifier for density model, model number in rows, parameters (y/n) in columns
  #visit_re is the visit random effect
  count.list <- matrix(NA, 4, 8)
  colnames(count.list) <- c("int", "year2013", "year2014", "year2015", "seasonSummer", "seasonSpring", "seasonAutumn", "visit_re")
  count.list[1, ] <- c(rep(1, 4), rep(0, 3), 1)     # year 
  count.list[2, ] <- c(1, 0, 0, 0, 1, 1, 1, 1)      # season
  count.list[3, ] <- c(rep(1, 8))                   # all variables
  count.list[4, ] <- c(1, rep(0, 6), 1)             # null model
  
  ################## picking the first model for detection function
  det.model[2] <- sample(1:nrow(det.list), size = 1)   # randomly choose a model to start with
  cur.dmod <- det.model[2]
  # holds the current det function parameters (vector \bmath{\theta}^t_m for model m)
  rj.cursigs <- det.param[2, ]
  
  # set proposal mean and sd
  det.prop.mean <- c(320, 3.6, rnorm(16, 0, 1))
  det.prop.sd <- c(20, 0.5, rep(0.1, 16))
  count.prop.mean <- c(1, rnorm(6, 0, 1), 0)
  count.prop.sd <- c(1, rep(0.1, 6), 1)
  
  ################## picking the first model for density model
  
  #count.model[2] <- sample(1:nrow(count.list), size = 1)   # randomly choose a model to start with
  count.model[2] <- 4
  cur.mod <- count.model[2]
  # holds the parameter values for the current density model (vector \bmath{\beta}^t for model m)
  rj.curparam <- count.param[2, ]
  
  #set prior min and max
  #currently uniform uninformative priors
  det.priors <- rbind(rep(-1, length(det.prop.mean)), rep(1, length(det.prop.mean)))
  count.priors <- rbind(rep(-1, length(count.prop.mean)), rep(1, length(count.prop.mean)))
  
  ############################### the priors ###########################################################
  
  
  ### detection function parameters
  # for scale intercept
  l.prior.sig <- function(sigm) {
    log.u.sig <- NULL
    for (k in 1:length(sigm)) {
      log.u.sig[k] <- log(dunif(sigm[k], 150, 500))                                
    }
    return(sum(log.u.sig))
  }
  
  # for shape
  l.prior.sha <- function(shap){
    log.u.sha <- NULL
    for (k in 1:length(shap)){
      log.u<-log(dunif(shap[k], 1, 8))
      if (is.infinite(log.u)) {
        log.u.sha[k]<- -100000
      } else {
        log.u.sha[k] <- log.u
      }
    }
    return(sum(log.u.sha))
  }
  
  # prior for random effect standard deviation (std.ran)
  l.prior.std.ran <- function (std.ran) {
    l.u.std.ran <- log(dunif(std.ran, 0, 2))
    if(is.infinite(l.u.std.ran)) {
      l.u.std.ran <- -100000
    }
    return(l.u.std.ran)
  }
  
  
  ### the priors for density model parameters
  
  l.prior <- function (x, min, max) {
    
    if (length(x) == 0) {
      return(0)
    }
    
    l.u_all <- NULL
    for (i in 1: length(x)) {
      l.u_all[i] <- log(dunif(x[i], min, max))
    }
    
    l.u  <- sum(l.u_all)
    if (is.infinite(l.u)) {
      l.u <- -100000
    }
    
    return(l.u)
    
  }
  
  
  #function to calculate likelihood for detection function for each observation
  calcFe <- function (x, sig.y, sig.s, sig.cc, sig.wc, sha2, sig1, sig.ss) {
    
    #calculate scale as a function of all parameters
    scale_param <- sig1 * exp(sig.y[years == x[8]] +
                                sig.s[x[2]] +
                                sig.ss[x[4]] +
                                sig.cc[x[5] + 1] +
                                sig.wc[x[6]])

    fe <- f.gamma.function(x[3], scale_param, sha2)
    

    return (fe)
    
  }
  
  #function to calculate the poisson likelihood of counts for each visit
  poissonLik <- function (x, int, b, yea, sea) {

    count_i <- counts[x["Visit"]] #count for visit i

    lambda <- exp(int + yea[years == x["Year"]] + sea[x["Season"]] + b[x["Visit"]])

    l.pois.y <- dpois(count_i, lambda)    # Poisson log-likelihood for each observation n_jpr in Y

    return(c(lambda, l.pois.y))

  }
  poissonLik <- cmpfun(poissonLik)
  
  
  ########################## the posterior conditional distribution functions ######################
  
  ### integrated likelihood L_{n,y}(\bmath{\beta},\bmath{\theta}} from eqn (1)  (LN: eqn. 2.1)
  # p combines all parameters \bmath{\theta}\bmath{\beta} and bj from eqns (1) - (7) (LN: 2.1-2.7)
  log.lik.fct <- function (p) {
    
    sig1 <- p[1]               # det fct: scale intercept
    sha2 <- p[2]               # det fct: shape
    sig.y <- c(p[3:4], 0)      # det fct: year 2013,2014,2015 coef
    sig.s <- c(p[5:6], 0)      # det fct: season Summer, Spring, Autumn coef
    sig.ss <- c(p[7:8], 0)     # det fct: sea_state 1, 2, 3 coef
    sig.cc <- c(p[9:16], 0)    # det fct: cloud_cover 0-8 coef
    sig.wc <- c(p[17:18], 0)   # det fct: water_clarity 1, 2, 3 coef
    
    int <- p[19]               # density: intercept
    yea <- p[20:22]            # density: year 2013, 2015, 2016
    sea <- p[23:25]            # density: season Summer, Spring, Autumn coef
    
    std.ran <- p[26]           # density: random effect standard deviation
    b <- p[27:length(p)]       # density: random effect coefficients

        #------------------------------------------------------------------------------
    # LN: Reasoning for 1 & 2 found in paragraph following eqn. 2.5
    #------------------------------------------------------------------------------
    
    
    # 3. calculate the f_e for each detection for det model likelihood component L_y(\bmath{\theta}) (eqn (3): exact distance data) (LN: eqn. 2.3)
    fe <- apply(cd, 1, calcFe, sig.y, sig.s, sig.cc, sig.wc, sha2, sig1, sig.ss)
    u <- integrate(f.gamma.function, truncate_left, truncate_right, sig1, sha2)$value/(truncate_right - truncate_left)    
    
    # 5. model L_n(\bmath{\beta}|\bmath{\theta}) from eqn (7)  (LN: eqn. 2.8)
    # for each visit 
    z <- rbinom(nz, 1, int) # latent indicator variables from data augmentation
    x <- runif(nz, truncate_left, truncate_right)
    p     <- f.gamma.function(x, sig1, sha2)
    mu   <- z * p

    #count likelihood
    count <- dbinom(round(mean(counts)), sum(z), mean(mu[mu != 0]))


    # pois_ll <- apply(visit_tab, 1, poissonLik, int, b, yea, sea)
    # lambda <- pois_ll[1, ]
    # l.pois.y <- pois_ll[2, ]
    l.b.norm <- log(dnorm(b, 0, std.ran))  # log of normal density for b_j
    
    #post <- log(prod((fe/u)/3)) + log(prod(l.pois.y)) + sum(l.b.norm)
    post <- log(prod((fe/u)/3)) + log(count) + sum(l.b.norm)
    
    if (is.infinite(post)) {
      post <- -100000
    }
    
    return(post)
  }
  
  
  #----------- DETECTION FUNCTIONS ---------#
  
  # using a half-normal detection function for line transects for f*y); g(y|\bmath{\theta}) from eqn (2) is given by (\pi(y) can be ommitted for line transects)
  f.hn.function <- function(dis, sigma) {
    f <- exp(-dis^2/(2*sigma^2))
    f
  }
  
  
  # using a hazard-rate function for line transects for f(y); g(y|\bmath{\theta}) from eqn (2) is given by (\pi(y) can be ommitted for line transects)
  f.haz.function <- function(dis, sigma, shape) {
    f <- 1-exp(-(dis/sigma)^(-shape))
    return(f)
  }
  
  #gamma detection function from the mrds package: keyfct.gamma
  f.gamma.function <- function (distance, scale, shape) {

    fr <- (1/gamma(shape)) * (((shape - 1)/exp(1))^(shape - 1))
    v1 <- distance/(scale * fr)
    return(v1^(shape-1)*exp(-v1)/(gamma(shape)*fr))

  }
  
  f.gamma.function <- cmpfun(f.gamma.function)
  
  
  #set up progress bar
  pb <- tkProgressBar(max=nt, title = paste("Chain", chain))
  progress <- function(n) {
    setTkProgressBar(pb, n)
  }
  
  #################################  the RJMCMC algorithm ######################################
  # nt is the number of iterations and is set above
  # row 1 is filled in with initial values for parameters and models
  for (i in 1:nt) {
    
    tryCatch({
      
      index <- i %% df_size + 2
      
      ##################### RJ step : sequential proposals to switch to another randomly selected model #####
      # all models are considered equally likely, i.e. P(m|m') = P(m'|m) for all m' and m
      
      
      ############## the detection function  ########################
      
      # the current model
      cur.dmod <- det.model[index - 1]
      curpa <- det.list[cur.dmod, ]
      
      #new_model <- sample(x = (1:nrow(det.list))[-cur.dmod], size = 1) #current model can't be chosen again
      new_model <- 8 #current model can't be chosen again
      newpa <- det.list[new_model, ]
      rj.newsigs <- rj.cursigs
      

      #fill in indeces for parameters that were not in the old model but are in the new one
      added_indeces <- curpa - newpa == -1
      rj.newsigs[added_indeces] <- rnorm(sum(added_indeces), det.prop.mean[added_indeces], det.prop.sd[added_indeces])

      #remove parameters that are in the old model but are not in the new one
      removed_indeces <- newpa - curpa == -1
      rj.newsigs[removed_indeces] <- 0

      num <- log.lik.fct(c(rj.newsigs, rj.curparam)) + l.prior(rj.newsigs[added_indeces], -1, 1) + sum(log(dnorm(rj.cursigs[removed_indeces], msyt.prop.mean[removed_indeces], msyt.prop.sd[removed_indeces])))  # the numerator of eqn (11)    (LN: Pretty sure this is A.4)
      den <- log.lik.fct(c(rj.cursigs, rj.curparam)) + sum(log(dnorm(rj.newsigs[added_indeces], msyt.prop.mean[added_indeces], msyt.prop.sd[added_indeces]))) + l.prior(rj.cursigs[removed_indeces], -1, 1) # the denominator of eqn (11)

      #check whether the new model is accepted
      A <- min(1, exp(num - den))                   # proposed move is accepted with probability A
      if (runif(1) <= A) {                           # if move is accepted change current values to new values
        rj.cursigs <- rj.newsigs
        cur.dmod <- new_model #change to new model if accepted

      } else {
        rj.newsigs <- rj.cursigs               # if move is rejected, reset everything to current
      }


      # record the model selection for the det fct in det.model for the i'th iteration
      det.model[index] <- cur.dmod
      
      
      #################### RJ step for density model ##################################
      rj.newparam <- rj.curparam
      
      # the current model:
      cur.mod <- count.model[index - 1]
      
      #choose a new model
      #new_model <- sample(x = (1:nrow(count.list))[-cur.mod], size = 1) #current model can't be chosen again
      new_model <- 4
      cur.par <- count.list[cur.mod, ]
      new.par <- count.list[new_model, ]
      
      
      #fill in indeces for parameters that were not in the old model but are in the new one
      added_indeces <- which(cur.par - new.par == -1)
      rj.newparam[added_indeces] <- rnorm(length(added_indeces), count.prop.mean[added_indeces], count.prop.sd[added_indeces])
      
      #remove parameters that are in the old model but are not in the new one
      removed_indeces <- which(new.par - cur.par == -1)
      rj.newparam[removed_indeces] <- 0
      
      num <- log.lik.fct(c(rj.cursigs, rj.newparam)) + l.prior(rj.newparam[added_indeces], -2, 2) + sum(log(dnorm(rj.newparam[removed_indeces], count.prop.mean[removed_indeces], count.prop.sd[removed_indeces])))
      den <- log.lik.fct(c(rj.cursigs, rj.curparam)) + sum(log(dnorm(rj.newparam[added_indeces], count.prop.mean[added_indeces], count.prop.sd[added_indeces]))) + l.prior(rj.newparam[removed_indeces], -2, 2)
      
      #check if the new model is accepted
      A <- min(1, exp(num - den))
      if (runif(1) <= A) {                             
        rj.curparam <- rj.newparam   #if accepted, update current model            
        cur.mod <- new_model
      } else {                             
        rj.newparam <- rj.curparam     #if rejected, reset parameters
      }
      
      # record which model we ended up with
      count.model[index] <- cur.mod
      
      ########################## Metropolis Hastings update ########################################################
      
      ########## updating the detection function parameters
      mh.newsigs <- rj.cursigs
      mh.cursigs <- rj.cursigs
      
      # for scale intercept
      mh.newsigs[1] <- mh.cursigs[1] + rnorm(1, 0, 2) 
      while (mh.newsigs[1] <=0) {    # prevents scale intercept to become < 0
        mh.newsigs[1] <- mh.cursigs[1] + rnorm(1, 0, 2)
      }
      num <- log.lik.fct(c(mh.newsigs, rj.curparam)) + l.prior.sig(mh.newsigs[1]) # the numerator of eqn (8)   (LN: This is is A.1)
      den <- log.lik.fct(c(mh.cursigs, rj.curparam)) + l.prior.sig(mh.cursigs[1]) # the denominator of eqn (8)
      A <- min(1, exp(num-den))
      if (runif(1) <= A & num != -1e5) {
        mh.cursigs <- mh.newsigs
      } else {
        mh.newsigs <- mh.cursigs
        
      }
      
      # for shape                  
      mh.newsigs[2] <- mh.cursigs[2] + rnorm(1, 0, 0.05) #tweaked for gamma to stop the shape parameter becoming < 1
      while (mh.newsigs[2] <=1) {
        mh.newsigs[2] <- mh.newsigs[2] + rnorm(1, 0, 0.01) 
      }
      num <- log.lik.fct(c(mh.newsigs, rj.curparam)) + l.prior.sha(mh.newsigs[2])
      den <- log.lik.fct(c(mh.cursigs, rj.curparam)) + l.prior.sha(mh.cursigs[2])
      A <- min(1, exp(num-den))
      if (runif(1) <= A & num != -1e5) {
        mh.cursigs <- mh.newsigs
      } else {
        mh.newsigs <- mh.cursigs
      }
      
      # for each coefficient in the scale parameter:
      indeces <- which(mh.cursigs[3:length(mh.cursigs)] != 0) + 2
      den_ll <- log.lik.fct(c(mh.cursigs, rj.curparam))
      
      for (ip in indeces) {
        
        mh.newsigs[ip] <- mh.cursigs[ip] + rnorm(1, 0, 0.12)
        num_ll <- log.lik.fct(c(mh.newsigs, rj.curparam))
        num <-  num_ll + l.prior(mh.newsigs[ip], -1, 1)
        den <-  den_ll + l.prior(mh.cursigs[ip], -1, 1)
        A <- min(1, exp(num - den))
        if (runif(1) <= A) {
          mh.cursigs <- mh.newsigs
          den_ll <- num_ll
        } else {
          mh.newsigs <- mh.cursigs
        } 
        
      }
      
      
      # fill in the new parameter values
      det.param[index, ] <- mh.cursigs
      rj.cursigs <- mh.cursigs
      
      
      
      ######### updating the density model parameters #############
      
      curparam <- rj.curparam
      newparam <- rj.curparam
      
      # the intercept
      newparam[1] <- curparam[1] + rnorm(1, 0, 0.01) 
      while (newparam[1] <=0 | newparam[1] >= 1) {
        newparam[1] <- curparam[1] + rnorm(1, 0, 0.01)
      }
      num <- log.lik.fct(c(rj.cursigs, newparam)) + l.prior(newparam[1], 0, 1)
      den <- log.lik.fct(c(rj.cursigs, curparam)) + l.prior(curparam[1], 0, 1) 
      A <- min(1, exp(num - den))
      if (runif(1) <= A & num != -1e5) {
        curparam[1] <- newparam[1]
      } else {
        newparam[1] <- curparam[1]
      }
      
      
     # loop through each parameter
      # indeces <- which(curparam[2:8] != 0) + 1
      # den_ll <- log.lik.fct(c(rj.cursigs, curparam)) 
      # for (m in indeces) {
      #   
      #   u <- rnorm(1, 0, 0.25)
      #   newparam[m] <- curparam[m] + u
      #   num_ll <- log.lik.fct(c(rj.cursigs, newparam))
      #   num <- num_ll + l.prior(newparam[m], -1, 1)
      #   den <- den_ll + l.prior(curparam[m], -1, 1)
      #   A <- min(1, exp(num - den))
      #   V <- runif(1)
      #   if (V <= A) {
      #     curparam[m] <- newparam[m]
      #     den_ll <- num_ll
      #   } else {
      #     newparam[m] <- curparam[m]
      #   }
        
        
      # }
      
      
      # the visit random effect standard deviation
      newparam[8] <- curparam[8] + max(rnorm(1, 0, 0.08), -newparam[8]) # cannot become 0 or less
      num <- log.lik.fct(c(rj.cursigs, newparam)) + l.prior.std.ran(newparam[8]) #changed log.ran.fct to log.lik.fct because not defined and probably the same
      den <- log.lik.fct(c(rj.cursigs, curparam)) + l.prior.std.ran(curparam[8]) #changed log.ran.fct to log.lik.fct because not defined and probably the same
      A <- min(1, exp(num-den))
      if (runif(1) <= A) {
        curparam[8] <- newparam[8]
      } else {
        newparam[8] <- curparam[8]
      } 
      
      # visit random effect coefficients
      
      den <- log.lik.fct(c(rj.cursigs, curparam))
      u <- rnorm(length(curparam), 0, 0.4)
      for (m in 9:length(curparam)) {
        newparam[m] <- curparam[m] + u[m]
        num <- log.lik.fct(c(rj.cursigs, newparam))
        A   <- min(1, exp(num-den))
        if (runif(1) <= A) {
          curparam[m] <- newparam[m]
          den <- num
        } else {
          newparam[m] <- curparam[m]
        }
      }
      
      #population estimate from data augmentation
      z <- rbinom(nz, 1, curparam[1]) # latent indicator variables from data augmentation
      N[index] <- sum(z)
      
      # saving the new parameter values of the density model in count.param
      count.param[index, ] <- curparam
      
      rj.curparam <- curparam
      rj.newparam <- curparam
      
      
      # saving the parameter matrices ever 1000 iterations
      if ((i + 1) %% df_size == 0) {
        
        print(paste("Chain", chain, "- Iteration", i + 1))
        
        progress(i) #iterate progress bar
        
        write.table(cbind(det.model, count.model), paste0("models", chain ,".csv"), append = TRUE, row.names = F, col.names = F, sep = ",")
        write.table(det.param, paste0("det.param", chain, ".csv"), append = TRUE, row.names = F, col.names = F, sep = ",")
        write.table(count.param, paste0("count.param", chain, ".csv"), append = TRUE, row.names = F, col.names = F, sep = ",")
        write.table(N, paste0("N", chain, ".csv"), append = TRUE, row.names = F, col.names = F, sep = ",")
        
        det.model[1] <- det.model[index]
        count.model[1] <- count.model[index]
        det.param[1, ] <- det.param[index, ]
        count.param[1, ] <- count.param[index, ]
        
      }
    }
    ,
    error = function(err) {
      print(paste("Chain", chain, ":Retrying run", i))
      i <- i - 1
    }
    )
  } # end of iteration 
  
}
