###########################################################################################################################

# Author: Cornelia Oedekoven 

# Centre for Research into Ecological and Environmental Modelling
# University of St Andrews, St Andrews, Scotland
# cornelia@mcs.st-and.ac.uk 

# R code for RJMCMC algorithm for analysing covey data (repeated point transects with exact distance data)

# Additional functions are given to adapt analysis to line transect and/or interval distance data 


###########################################################################################################################

# Data format required:

# for the detection function model L_y(\bmath{\theta}) (eqn (2.3))
# matrix covey.d: a n*6 matrix  (n = total number of detections, 6 columns: id, distance, site, year, type, state)
# LN: id and distance will be retained, the rest should be relevent covariates to the shark transects


if (Sys.info()[4] == "SCI-6246") {
  setwd(dir = "C:/Users/43439535/Documents/Lisa/phd/aerial survey/data")
  source_location <- "~/Lisa/phd/Mixed models/R code/R-functions-southern-ocean/"
} else {
  setwd(dir = "C:/Users/Lisa/Documents/phd/aerial survey/data")
  source_location <- "~/phd/southern ocean/Mixed models/R code/R-functions-southern-ocean/"
}

#covey.d<-covey2[which(is.na(covey2$Distance)==F),]
dat <- read.csv("aerial_survey_summary_r.csv", header = T) #lisa's sighting data
covey.d <- dat[dat$Species == "BOT" & !is.na(dat$Dist..from.transect) & dat$Secondary != "Y", c("Date", "Season", "Dist..from.transect", "Beaufort.Sea.State", "Cloud.cover", "Water.clarity")]
covey.d$Date <- as.character(covey.d$Date)
covey.d$Year <- as.numeric(substr(covey.d$Date, nchar(covey.d$Date) - 3, nchar(covey.d$Date))) #extract year from date
names(covey.d) <- c("Date", "Season", "Distance", "Sea_state", "Cloud_cover", "Water_clarity", "Year")
covey.d$Id <- 1:nrow(covey.d)
covey.d$Visit <- as.numeric(as.factor(covey.d$Date)) #visit is transect
covey.d <- covey.d[covey.d$Distance != 0, ] #remove 0 distances because they are errors

length.d <- length(covey.d$Distance)

#-----------------------------------------------------------------
# LN: This section is mostly specific to the case study in the ms
#     We will only have one site and no points
#     Need to identify any potential covariates
#-----------------------------------------------------------------

# for the density model L_n(\bmath{\beta}|\bmath{\theta}) (eqn (6)) (LN: Start of Sec. 2.2)
# matrix n.jpr: a j*max(Tj) matrix with n_jpr (counts at site j, point p, visit r (Tj is a vector of numbers of measurements in the jth group)
# matrices with covariates: j*max(Tj) matrix with observed covariate values at site j, point p, visit r: Year, Type, Day (Julian day) and State

# number of observations for each site
#LN: We only have one site, so this may not be needed

#site<-sort(unique(count.data$Site))
#j<-length(site)
years          <- sort(unique(covey.d$Year))
seasons        <- levels(covey.d$Season)
sea_states     <- sort(unique(covey.d$Sea_state))
cloud_covers   <- sort(unique(covey.d$Cloud_cover))
water_claritys <- sort(unique(covey.d$Water_clarity))
visits         <- sort(unique(covey.d$Visit))

j  <- length(visits)
Y <- table(covey.d$Visit) #count per visit


# matrix that will hold the counts n.jpr for each visit (1 row per site)
#Y <- matrix(NA, j, max(Tj))

# matrix that will hold the values for Year (1 for 2006, 2 for 2007, 3 for 2008): factor covariate with 2 levels
#Year<-matrix(NA,j,max(Tj))

# matrix that will hold a 0 or 1 depending on whether Type = CONTROL or TREAT, respectively
#Type<-matrix(NA,j,max(Tj))

# matrix that will hold the values for Julian day
#Day<-matrix(NA,j,max(Tj))

# matrix that will hold the values for State
#State<-matrix(NA,j,max(Tj))

# filling in the above values from
# for (i in 1:j){
#   
#   x<-which(count.data$Site==site[i])
#   l<-length(count.data$Count[x])
#   y<-count.data$Count[x]
#   for (k in 1:l){
#     
#     Y[i,k]<-y[k]
#     Year[i,k]<-ifelse(count.data$Year[x[k]]==2006,1,{ifelse(count.data$Year[x[k]]==2007,2,3)})
#     Type[i,k]<-ifelse(count.data$Type[x[k]]=="TREAT",1,0)
#     Day[i,k]<-count.data$jd[x[k]]
#     State[i,k]<-which(States==count.data$State[x[k]])
#   }
# }
# sum(Y[which(is.na(Y)==F)])
#[1] 2545

########################################  set initial values ##############

# global hazard-rate model with scale and shape for L_y(\bmath{\theta}) (eqn 2.3)
scale0 <- 130
shape0 <- 2.5

#  count model parameters: intercept and random effect standard deviation for L_n(\bmath{\beta}|\bmath{\theta}) (eqn (6))
int0 <- -13
std.ran0 <- 1

# the random effect coefficients b_j
# LN: We only have one site, so can't really use this
#could try using visit as random effect
#b0<-rnorm(j,0,std.ran0)


#########################################################################
# setting up the matrices that will contain the paramter values;

# number of iterations
#nt <- 100000
nt <- 10 #fewer iterations for testing

#-------------------------------------------------------------------------
# LN: Indexing here is specific to the case study
#     Will need to adapt it to match the shark data
#     Remember that parameters include standard deviation estimates
#-------------------------------------------------------------------------

# holds the values for detection function parameters for each iteration
# 15 colums due to 15 parameters in full model: hazard-rate det fct with covariates: year, type and state
# we have 23 parameters: scale, shape, year(3), season(3), sea_state(3), cloud(9), water_clarity(3)
det.param <- matrix(NA, nt+1, 18)

# for an intercept only model:
det.param[1, ] <- c(scale0, shape0, rep(0, 16))

# holds the model id number for detection function for each iteration, 
det.model <- matrix(NA, nt+1, 1)            # refers to det.list

# the matrix that will keep the parameter values for the density model
count.param <- matrix(NA, nt+1, 22)    

# filling in the initial values
count.param[1, 1:22] <- c(int0, rep(0, 21))

# holds the model id number for density model for each iteration
count.model <- matrix(NA, nt+1, 1)          # refers to count.list

############# proposal distributions
# proposal distributions for detection function parameters:

# 1. for main analysis and prior sensitivity analysis
det.prop.mean <- c(138.60, 3.00, rnorm(16, 0, 1))
det.prop.sd <- c(1.41, 0.84, rep(0.1, 16))

# proposal distribution for the fixed effect density model parameters
# 1. for covey data
count.prop.mean <- c(-13.06, 0, rnorm(21, 0, 1))
count.prop.sd <- c(0.30, 0, rep(0.1, 21))

msyt.prop.mean <- c(1, rep(0.5, 21))
msyt.prop.sd <- rep(0.5, 22)

################## picking the first model for detection function
det.model[1, 1] <- 5   # global hazard-rate: \bmath{\theta} = \{\sigma,\tau\}
cur.dmod <- det.model[1, 1]
# holds the current det function parameters (vector \bmath{\theta}^t_m for model m)
rj.cursigs <- det.param[1, ]


################## picking the first model for density model
# there are 16 models (all of them include random effect for Pair2):
# to pick the first model:
count.model[1] <- 1          # intercept and random effect only model \bmath{\beta}=\{\beta_0,\sigma_b^2\}
cur.mod <- 1
# holds the parameter values for the current density model (vector \bmath{\beta}^t for model m)
rj.curparam <- count.param[1, 1:22]

#-------------------------------------------------------------------
# LN: This will change since we will define our own models
#-------------------------------------------------------------------
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
count.list <- matrix(NA, 8, 22)
colnames(count.list)<-c("int", "year2013", "year2014", "year2015", "seasonSummer", "seasonSpring", "seasonAutumn", paste0("ss", 1:3),  paste0("cc", 0:8), paste0("wc", 1:3))
count.list[1, ] <- c(rep(1, 4), rep(0, 18))                # year 
count.list[2, ] <- c(1, 0, 0, 0, 1, 1, 1, rep(0, 15))      # season
count.list[3, ] <- c(rep(1, 7), rep(0, 15))                # year and season
count.list[4, ] <- c(1, rep(0, 6), 1, 1, 1, rep(0, 12))    # sea state
count.list[5, ] <- c(1,rep(0,21))                          # global hazard rate model
count.list[6, ] <- c(1,rep(1,21))                          # all variables
count.list[7, ] <- c(1, rep(0, 18), 1, 1, 1)               # water clarity
count.list[8, ] <- c(1, rep(0, 9), rep(1, 9), 0, 0, 0)     # cloud cover

################################# Likelihood functions ################################################

# hazard-rate function for point transects: \pi(y) * g(y|\bmath{\theta}) from eqn (2) (LN: Paragraph following eqn. 2.2)
f.haz.function<-function(dis, sigma, shape) {
  f <- 2*pi*dis*(1-exp(-(dis/sigma)^(-shape)))
  return(f)
}

############################### the priors ###########################################################
#-----------------------------------------------------------------
# LN: Can use the some of the existing functions
#     Need to make sure limits are reasonable for the shark data
#     Will need to make sure the shark covariates are represented
#-----------------------------------------------------------------


### detection function parameters
# for scale intercept
l.prior.sig <- function(sigm) {
  log.u.sig<-array(NA,length(sigm))
  for (k in 1:length(sigm)) {
    log.u.sig[k] <- log(dunif(sigm[k], 1, 100000))                                
  }
  return(sum(log.u.sig))
}

# for shape
l.prior.sha<-function(shap){
  log.u.sha<-array(NA,length(shap))
  for (k in 1:length(shap)){
    log.u<-log(dunif(shap[k], 1, 20))
    ifelse(abs(log.u) == Inf,log.u.sha[k]<- -100000, log.u.sha[k] <- log.u)}
  return(sum(log.u.sha))
}

#for scale coefficients
l.prior.coef<-function(coefsig){
  lcs<-length(coefsig)
  log.u.coefsig<-array(NA,lcs)
  for (k in 1:lcs){
    log.u<-log(dunif(coefsig[k],-2.5,2.5))
    ifelse(abs(log.u)==Inf,log.u.coefsig[k]<- -100000,log.u.coefsig[k]<-log.u)}
  return(sum(log.u.coefsig))}

# prior for random effect standard deviation (std.ran)
l.prior.std.ran <- function (std.ran) {
  l.u.std.ran <- log(dunif(std.ran, 0, 2))
  return(l.u.std.ran)
}


### the priors for density model parameters

l.prior <- function (x, min, max) {
  
  #x = quantiles for dunif. I per level of factor
  
  l.u_all <- NULL
  for (i in 1: length(x)) {
    l.u_all[i] <- log(dunif(x[i], min, max))
  }
  
  l.u  <- sum(l.u_all)
  if (abs(l.u) == Inf) {
    l.u <- -100000
  }
  
  return(l.u)
  
}


########################## the posterior conditional distribution functions ######################

### integrated likelihood L_{n,y}(\bmath{\beta},\bmath{\theta}} from eqn (1)  (LN: eqn. 2.1)
# p combines all parameters \bmath{\theta}\bmath{\beta} and bj from eqns (1) - (7) (LN: 2.1-2.7)
log.lik.fct <- function (p) {
  
  sig1 <- p[1]               # det fct: scale intercept
  sha2 <- p[2]               # det fct: shape
  sig.y <- c(p[3:4], 0)       # det fct: year 2013,2014,2015 coef
  sig.s <- c(p[5:6], 0)      # det fct: season Summer, Spring, Autumn coef
  sig.ss <- c(p[7:8], 0)      # det fct: sea_state 1, 2, 3 coef
  sig.cc <- c(p[9:16], 0)     # det fct: cloud_cover 0-8 coef
  sig.wc <- c(p[17:18], 0)    # det fct: water_clarity 1, 2, 3 coef
  
  int <- p[19]               # density intercept
  yea <- p[20:22]            # density year 2013, 2015, 2016
  sea <- p[23:25]            # density season Summer, Spring, Autumn coef
  ss  <- p[26:28]            # density sea_state 1, 2, 3 coef
  cc  <- p[29:37]            # density cloud_cover 0-8 coef
  wc  <- p[38:40]            # density water_clarity 1, 2, 3 coef
  
  #------------------------------------------------------------------------------
  # LN: Reasoning for 1 & 2 found in paragraph following eqn. 2.5
  #     Dimensions will need to be made specific to the shark data at each step
  #------------------------------------------------------------------------------
  # 1. calculate the different scale parameters as function of parameters
  
  combns <- expand.grid("year"= c("2013", "2014", "2015"), "season"= c("summer", "spring", "autumn"), "ss"= 1:3, 
                        "cc"= 0:8, "wc"= 1:3)  
  
  sig.msyt <- rep(NA, nrow(combns))      # 11 states (rows), 3 years * 2 type levels (CONTROL,TREAT) (columns)
  efa      <- rep(NA, nrow(combns))
  
  
  # for (strat in 1:11){        
  #   for (ty in 4:6){
  #     sig.msyt[strat,ty]<-sig1*exp(sig.st[strat]+ sig.y[ty-3])}}
  # 
  # for (strat in 1:11){
  #   for (ty in 1:3){
  #     sig.msyt[strat,ty]<-sig.msyt[strat,ty+3]*exp(sig.t)}}
  
  
  sig.msyt <- sig1 * exp(sig.y[as.numeric(combns[, "year"])] + sig.s[as.numeric(combns[, "season"])] + 
                           sig.ss[as.numeric(combns[, "ss"])] + sig.cc[as.numeric(combns[, "cc"]) + 1] +
                           sig.wc[as.numeric(combns[, "wc"])])
  
  
  
  # 2. calculate the different effective areas as a function of covariates (using the scales from sig.msyt)
  # for (strat in 1:11){
  #   for (ty in 1:6){
  #     efa[strat,ty]<-integrate(f.haz.function,0,500,sig.msyt[strat,ty],sha2)$value}}
  
  for (i in 1:length(sig.msyt)) {
     tryCatch (
     efa[i] <- integrate(f.haz.function, 0, 500, sig.msyt[i], sha2)$value,
    error = function(err) {
      print(p)})
  }
  
  # 3. calculate the f_e for each detection for det model likelihood component L_y(\bmath{\theta}) (eqn (3): exact distance data) (LN: eqn. 2.3)
  fe <- array(NA, length.d)
  offset <- array(NA, length.d)
  for (i in 1:length.d){
    
    combn_row <- which(combns$year == covey.d$Year[i] & combns$season == covey.d$Season[i] &
                         combns$ss == covey.d$Sea_state[i] & combns$cc == covey.d$Cloud_cover[i] &
                         combns$wc == covey.d$Water_clarity[i])
    
    norm.const <- efa[combn_row]    # normalising constant from denominator in eqn (2) equals the effective area for a given scale and shape parameter
    fe[i] <- log(f.haz.function(covey.d$Distance[i], sig.msyt[combn_row], sha2)/norm.const)
    
    offset[i] <- norm.const
    
  }
  
  # 5. model L_n(\bmath{\beta}|\bmath{\theta}) from eqn (7)  (LN: eqn. 2.8)
  l.pois.y <- NULL  # matrix that will hold the Poisson likelihood for each observation n_jpr
  lambda   <- NULL     # matrix for storing the lambda_jpr from eqn (6)
  # for each visit 
  for (i in visits) {
    
    count_i <- sum(covey.d$Visit == i) #count for visit i
    
    if (count_i == 0) {
      lambda[i] <- exp(int) #if 0 count on visit i, use only intercept
    } else {
      lambda[i] <- exp(int + yea[which(years == unique(covey.d$Year[covey.d$Visit == i]))] + sea[which(seasons == unique(covey.d$Season[covey.d$Visit == i]))] + log(mean(offset[covey.d$Visit == i])))
    }
    l.pois.y[i] <- log(dpois(count_i, lambda[i]))    # Poisson log-likelihood for each observation n_jpr in Y
  }
  
  post <- sum(fe) + sum(l.pois.y[!is.na(l.pois.y)])
  return(post)
}


################################# other function you will need
# this function matches a string of numbers against rows a matrix
# and will return the row number of matrix that matches the string exactly
match.function <- function (xmod,model.matrix) {
  is.match.or.not <- array(NA, length(model.matrix[, 1]))
  for (i in 1:length(model.matrix[, 1])) { 
    is.match.or.not[i] <- sum(xmod == model.matrix[i, ])}
  result <- which(is.match.or.not == length(model.matrix[1, ]))
  return(result)
}

#-------------------------------------------------------------------------------
# LN: We want to be using the functions for line transects
#-------------------------------------------------------------------------------

#need to change to gamma detection function at some stage

###### Alterations to this algorithm: 
# dis = distance (perpendicular for lines and radial for points, sigma = scale parameters, shape = shape parameter

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


#################################  the RJMCMC algorithm ######################################
# nt is the number of iterations and is set above
# row 1 is filled in with initial values for parameters and models
for (i in 2:nt) {
  print(i)
  
  ##################### RJ step : sequential proposals to add or delete covariates depending on whether they are in the model or not #####
  # all models are considered equally likely, i.e. P(m|m') = P(m'|m) for all m' and m
  
  ############## the detection function  ########################
  #--------------------------------------------------
  # LN: More adjusting for the relevent shark models
  #--------------------------------------------------
  
  # the current model
  cur.dmod <- det.model[i - 1]
  curpa <- det.list[cur.dmod, ]
  # setting the parameters for the new model equal to the current model
  newpa <- curpa
  rj.newsigs <- rj.cursigs
  
  # going through the list of coefficients to check whether to add or remove one
  
  for (param in c("year", "season", "ss", "cc", "wc")) {
    
    indeces <- grep(param, names(curpa))
    
    
    if (sum(curpa[indeces]) == 0) {        # if param is not currently in the model, propose to add it
      newpa[indeces] <- 1
      rj.newsigs[indeces] <- rnorm(length(indeces), det.prop.mean[indeces], det.prop.sd[indeces])   # draw random samples from proposal distributions
      num <- log.lik.fct(c(rj.newsigs, rj.curparam)) + l.prior.coef(rj.newsigs[indeces]) # the numerator of eqn (11)    (LN: Pretty sure this is A.4)
      den <- log.lik.fct(c(rj.cursigs, rj.curparam)) + sum(log(dnorm(rj.newsigs[indeces], msyt.prop.mean[indeces], msyt.prop.sd[indeces]))) # the denominator of eqn (11)
    } else {
      newpa[indeces] <- 0               # if param is in the current model, propose to delete it
      rj.newsigs[indeces] <- 0
      num <- log.lik.fct(c(rj.newsigs, rj.curparam)) + sum(log(dnorm(rj.cursigs[indeces], msyt.prop.mean[indeces], msyt.prop.sd[indeces]))) # the numerator of eqn (11)
      den <- log.lik.fct(c(rj.cursigs, rj.curparam)) + l.prior.coef(rj.cursigs[indeces]) # the denominator of eqn (11)
    }
    
    #check whether the new model is accepted
    A <- min(1, exp(num - den))                   # proposed move is accepted with probability A
    V <- runif(1)
    if (V <= A) {                           # if move is accepted change current values to new values
      rj.cursigs <- rj.newsigs               
      curpa <- newpa
    } else {                             
      rj.newsigs<-rj.cursigs               # if move is rejected, reset everything to current
      newpa<-curpa
    }
    
  }
  
  # which model did we end up with 
  cur.dmod <- match.function(curpa, det.list)
  
  #if model doesn't exist in our subset of chosen models, go back to the previous model
  #probably not the right thing to do in this situation
  if (length(cur.dmod) == 0) {
    cur.dmod <- det.model[i - 1]
  }
  
  # record the model selection for the det fct in det.model for the i'th iteration
  det.model[i] <- cur.dmod
  
  
  #################### RJ step for density model ##################################
  rj.newparam <- rj.curparam
  # the current model:
  cur.mod <- count.model[i - 1]
  # the parameters in the current model
  cur.par <- count.list[cur.mod, ]
  new.par <- cur.par
  
  # for each parameter, check if it is in the current model 
  
  for (param in c("year", "season", "ss", "cc", "wc")) {
    
    indeces <- grep(param, names(cur.par))
    
    if (sum(cur.par[indeces]) == 0) {                  # if param is not in current model, propose to add it
      new.par[indeces] <- 1
      
      for (f in indeces) {
        rj.newparam[f] <- rnorm(1,count.prop.mean[f], count.prop.sd[f])
      }
      
      rj.newparam[indeces] <- rnorm(1, count.prop.mean[indeces], count.prop.sd[indeces])
      num <- log.lik.fct(c(rj.cursigs, rj.newparam)) + l.prior(rj.newparam[indeces], -1, 1) 
      den <- log.lik.fct(c(rj.cursigs, rj.curparam)) + sum(log(dnorm(rj.newparam[indeces], count.prop.mean[indeces], count.prop.sd[indeces])))
      
      #check if the new model is accepted
      A <- min(1, exp(num - den))
      V <- runif(1)
      if (V <= A) {                             
        rj.curparam <- rj.newparam   #if accepted, update current model            
        cur.par <- new.par
      } else {                             
        rj.newparam <- rj.curparam     #if rejected, reset parameters
        new.par <- cur.par
      }
    } else {
      rj.newparam[indeces] <- 0         # if param is in the current model, propose to delete it
      new.par[indeces] <- 0
      num <- log.lik.fct(c(rj.cursigs,rj.newparam))  + sum(log(dnorm(rj.curparam[indeces],count.prop.mean[indeces],count.prop.sd[indeces])))
      den <- log.lik.fct(c(rj.cursigs,rj.curparam)) + l.prior(rj.curparam[indeces], -1, 1) 
      A <- min(1, exp(num - den))
      V <- runif(1)
      if ( V <= A) {                             
        rj.curparam <- rj.newparam               
        cur.par <- new.par
      } else {                             
        rj.newparam <- rj.curparam
        new.par <- cur.par
      }
    }  
  }
  
  # which model did we end up with:
  cur.mod <- match.function(cur.par, count.list)
  
  #if chosen model is not in the our list of possible models, go back to previous model
  if (length(cur.mod) == 0) {
    cur.mod <- count.model[i - 1]
  }
  
  count.model[i] <- cur.mod
  
  ########################## Metropolis Hastings update ########################################################
  
  ########## updating the detection function parameters
  mh.newsigs <- rj.cursigs
  mh.cursigs <- rj.cursigs
  
  # for scale intercept
  u <- rnorm(1, 0, 3.5)                
  if ((mh.cursigs[1] + u) > 1) {    
    # prevents scale intercept to become < 0
    mh.newsigs[1] <- mh.cursigs[1] + u
    # the numerator of eqn (8)   (LN: This is is A.1)
    num <- log.lik.fct(c(mh.newsigs, rj.curparam)) + l.prior.sig(mh.newsigs[1])
    # the denominator of eqn (8)
    den <- log.lik.fct(c(mh.cursigs, rj.curparam)) + l.prior.sig(mh.cursigs[1])
    A <- min(1, exp(num-den))
    V <- runif(1)
    ifelse(V <= A, mh.cursigs <- mh.newsigs, mh.newsigs <- mh.cursigs)    
  }
  
  # for shape                  
  u <- rnorm(1, 0, 0.2)
  mh.newsigs[2] <- mh.cursigs[2] + u
  num <- log.lik.fct(c(mh.newsigs, rj.curparam)) + l.prior.sha(mh.newsigs[2])
  den <- log.lik.fct(c(mh.cursigs, rj.curparam)) + l.prior.sha(mh.cursigs[2])
  A <- min(1, exp(num-den))
  V <- runif(1)
  ifelse(V <= A, mh.cursigs <- mh.newsigs, mh.newsigs <- mh.cursigs)
  
  # for each coefficient in the scale parameter:
  
  for (param in c("year", "season", "ss", "cc", "wc")) {
    
    indeces <- grep(param, names(mh.cursigs))
    
    if (sum(mh.cursigs[indeces]) != 0) { #only update parameters in the current model
      for (ip in indeces) {
        #u<-rnorm(1,0,0.01)          # acc prob = 80-95%
        u <- rnorm(1,0,0.12)
        mh.newsigs[ip]<-mh.cursigs[ip] + u
        num <- log.lik.fct(c(mh.newsigs, rj.curparam)) + l.prior.coef(mh.newsigs[ip])
        den <- log.lik.fct(c(mh.cursigs, rj.curparam)) + l.prior.coef(mh.cursigs[ip])
        A <- min(1, exp(num - den))
        V <- runif(1)
        ifelse(V <= A, mh.cursigs <- mh.newsigs, mh.newsigs <- mh.cursigs)
      }
    }
  }
  
  # fill in the new parameter values
  det.param[i, ] <- mh.cursigs
  rj.cursigs <- mh.cursigs
  
  
  
  ######### updating the density model parameters #############
  
  curparam <- rj.curparam
  newparam <- rj.curparam
  
  # the intercept
  u <- rnorm(1,0,0.08)                        
  newparam[1] <- curparam[1] + u
  num <- log.lik.fct(newparam) + l.prior(newparam[1], -20, 7)
  den <- log.lik.fct(curparam) + l.prior(curparam[1], -20, 7)
  A <- min(1, exp(num - den))
  V <- runif(1)
  ifelse(V <= A, curparam[1] <- newparam[1], newparam[1] <- curparam[1])
  
  # loop through each parameter
  
  for (param in c("year", "season", "ss", "cc", "wc")) {
    
    indeces <- grep(param, names(cur.par))
    
    if (sum(curparam[indeces]) == 0){ #only update parameters in the current model
      for (m in indeces){
        u<-rnorm(1, 0, 0.25)
        newparam[m] <- curparam[m] + u
        num <- log.lik.fct(c(rj.cursigs,newparam)) + l.prior(newparam[m], -1, 1)
        den <- log.lik.fct(c(rj.cursigs,curparam)) + l.prior(curparam[m], -1, 1)
        A <- min(1, exp(num - den))
        V <- runif(1)
        ifelse(V <= A, curparam[m] <- newparam[m], newparam[m] <- curparam[m])
      } 
    }
  }
  
  # the random effect standard deviation
  #change 18 to correct index for random effect standard deviation
  sd_index <- length(newparam)
  u <- max(rnorm(1, 0, 0.08), -newparam[sd_index])    # cannot become 0 or less
  newparam[sd_index] <- curparam[sd_index]+u
  num <- log.lik.fct(c(rj.cursigs, newparam)) + l.prior.std.ran(newparam[sd_index]) #changed log.ran.fct to log.lik.fct because not defined and probably the same
  den <- log.lik.fct(c(rj.cursigs, curparam)) + l.prior.std.ran(curparam[sd_index]) #changed log.ran.fct to log.lik.fct because not defined and probably the same
  A <- min(1, exp(num-den))
  V <- runif(1)
  ifelse(V <= A, curparam[sd_index] <- newparam[sd_index], newparam[sd_index] <- curparam[sd_index])
  
  #commented out for now but will add random effect for visit later
  # the random effects coefficients
  # for (m in 19:(j+18)){
  #   u<-rnorm(1,0,0.4)
  #   newparam[m]<-curparam[m]+u
  #   num<-log.lik.fct(c(rj.cursigs,newparam))
  #   den<-log.lik.fct(c(rj.cursigs,curparam))
  #   A<-min(1,exp(num-den))
  #   V<-runif(1)
  #   ifelse(V<=A,curparam[m]<-newparam[m],newparam[m]<-curparam[m])
  # }
  
  
  # saving the new parameter values of the density model in count.param
  count.param[i, ] <- curparam
  
  rj.curparam <- curparam
  rj.newparam <- curparam
  
  # saving the parameter matrices ever 1000 iterations
  if (i %% 1000 == 0) {
    save(det.model, file = 'det.model.RData')
    save(count.model, file = 'count.model.RData')
    save(det.param, file = 'msyt.param.RData')
    save(count.param, file = 'count.param.RData')
  }
  
} # end of iteration


