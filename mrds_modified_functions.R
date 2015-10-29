dht <- function(model,region.table,sample.table, obs.table=NULL, subset=NULL,
                se=TRUE, bootstrap=FALSE, options=list()){
  # Functions Used:  assign.default.values, create.varstructure,
  #                  covered.region.dht, survey.region.dht, dht.se, varn,
  #                  covn(in varn.R), solvecov (in coef.ds.R).
  
  tables.dht <- function(group){
    # Internal function to create summary tables for clusters (group=TRUE) or
    # individuals (group=FALSE).
    options$group <- group
    
    # Compute covered region abundances by sample depending on value of group
    Nhat.by.sample <- covered.region.dht(obs, samples, group)
    
    # Mod 18-Aug-05 jll; added computation of avergage detection probability
    # which is simply n/Nhat in the covered region
    average.p <- nrow(obs)/sum(Nhat.by.sample$Nhat)
    
    # Scale up abundances to survey region
    # jll 19-Jan-05 - sort Nhat.by.sample by Region.Label and Sample.Label
    width <- model$meta.data$width * options$convert.units
    Nhat.by.sample <- survey.region.dht(Nhat.by.sample, samples,width,point)
    Nhat.by.sample <- Nhat.by.sample[order(Nhat.by.sample$Region.Label,
                                           Nhat.by.sample$Sample.Label),]
    if(point){
      s.area <- Nhat.by.sample$Effort.x*pi*width^2
    }else{
      s.area <- Nhat.by.sample$Effort.x*width ####removed *2 to make one sided####
    }
    
    bysample.table <- with(Nhat.by.sample,
                           data.frame(Region      = Region.Label,
                                      Sample      = Sample.Label,
                                      Effort      = Effort.x,
                                      Sample.Area = s.area,
                                      Area        = Area,
                                      n           = n,
                                      Nhat        = Nhat,
                                      Nchat       = Nhat*CoveredArea/Area))
    
    bysample.table$Dhat <- bysample.table$Nchat/bysample.table$Sample.Area
    Nhat.by.region <- by(Nhat.by.sample$Nhat, Nhat.by.sample$Region.Label,sum)
    
    # Create estimate table
    numRegions <- length(unique(samples$Region.Label))
    if(numRegions > 1){
      estimate.table <- data.frame(
        Label = c(levels(unique(samples$Region.Label)),"Total"),
        Estimate = rep(0, numRegions + 1),
        se = rep(NA,numRegions + 1),
        cv = rep(NA, numRegions + 1),
        lcl = rep(NA,numRegions + 1),
        ucl = rep(NA, numRegions + 1))
    }else{
      estimate.table = data.frame(Label    = c("Total"),
                                  Estimate = rep(0,1),
                                  se       = rep(NA, 1),
                                  cv       = rep(NA, 1),
                                  lcl      = rep(NA, 1),
                                  ucl      = rep(NA, 1))
    }
    
    if(numRegions > 1){
      estimate.table$Estimate <- c(Nhat.by.region, sum(Nhat.by.region))
    }else{
      estimate.table$Estimate <- Nhat.by.region
    }
    
    # Create summary table
    summary.table <- Nhat.by.sample[, c("Region.Label","Area",
                                        "CoveredArea","Effort.y")]
    summary.table <- unique(summary.table)
    var.er <- sapply(split(Nhat.by.sample,Nhat.by.sample$Region.Label),
                     function(x)varn(x$Effort.x,x$n,type=options$ervar))
    
    if(numRegions > 1){
      var.er <- c(var.er,varn(Nhat.by.sample$Effort.x,
                              Nhat.by.sample$n,type=options$ervar))
    }
    
    #  jll 11_11_04; change to set missing values for nobs to 0
    #   - regions with no sightings
    nobs <- as.vector(by(bysample.table$n, bysample.table$Region, sum))
    nobs[is.na(nobs)] <- 0
    summary.table$n <- nobs
    if(group){
      summary.table$k <- tapply(Nhat.by.sample$Sample.Label,
                                Nhat.by.sample$Region.Label,length)
      colnames(summary.table) <- c("Region", "Area", "CoveredArea",
                                   "Effort", "n","k")
    }else{
      colnames(summary.table) = c("Region", "Area", "CoveredArea", "Effort", "n")
    }
    
    if(numRegions > 1){
      summary.table <- data.frame(Region=c(levels(summary.table$Region),"Total"),
                                  rbind(summary.table[, -1],
                                        apply(summary.table[,-1], 2, sum)))
    }
    
    summary.table$ER <- summary.table$n/summary.table$Effort
    summary.table$se.ER <- sqrt(var.er)
    summary.table$cv.ER <- summary.table$se.ER/summary.table$ER
    summary.table$cv.ER[summary.table$ER==0] <- 0
    
    # set missing values to 0
    summary.table$ER[is.nan(summary.table$ER)] <- 0
    summary.table$se.ER[is.nan(summary.table$se.ER)] <- 0
    summary.table$cv.ER[is.nan(summary.table$cv.ER)] <- 0
    
    # If summary of individuals for a clustered popn, add mean
    #  group size and its std error
    if(!group){
      mean.clustersize <- tapply(obs$size,obs$Region.Label,mean)
      se.clustersize <- sqrt(tapply(obs$size,obs$Region.Label,var)/
                               tapply(obs$size,obs$Region.Label,length))
      cs <- data.frame(Region=names(mean.clustersize),
                       mean.size=as.vector(mean.clustersize),
                       se.mean=as.vector(se.clustersize))
      
      summary.table <- merge(summary.table, cs, by.x = "Region",
                             all=TRUE,sort=FALSE)
      if(numRegions > 1){
        summary.table$mean.size[numRegions+1] <- mean(obs$size)
        summary.table$se.mean[numRegions+1] <-sqrt(var(obs$size)/length(obs$size))
      }
      #  29/05/12 lhm - moved to set missing values to 0
      summary.table$mean.size[is.na(summary.table$mean.size)] <- 0
      summary.table$se.mean[is.na(summary.table$se.mean)] <- 0
    }
    
    rownames(summary.table) <- 1:dim(summary.table)[1]
    
    # If a std error has been requested call dht.se
    if(se){
      result <- dht.se(model, summary.table, samples, obs, options, numRegions,
                       estimate.table, Nhat.by.sample)
      estimate.table <- result$estimate.table
    }
    
    # Create estimate table for D from same table for N
    D.estimate.table <- estimate.table
    if(numRegions > 1){
      D.estimate.table$Estimate <-D.estimate.table$Estimate/c(region.table$Area,
                                                              sum(region.table$Area))
      D.estimate.table$se <- D.estimate.table$se/c(region.table$Area,
                                                   sum(region.table$Area))
      D.estimate.table$lcl <- D.estimate.table$lcl/c(region.table$Area,
                                                     sum(region.table$Area))
      D.estimate.table$ucl <- D.estimate.table$ucl/c(region.table$Area,
                                                     sum(region.table$Area))
    }else{
      D.estimate.table$Estimate <- D.estimate.table$Estimate/region.table$Area
      D.estimate.table$se <- D.estimate.table$se/region.table$Area
      D.estimate.table$lcl <- D.estimate.table$lcl/region.table$Area
      D.estimate.table$ucl <- D.estimate.table$ucl/region.table$Area
    }
    
    # set missing values to 0
    D.estimate.table$Estimate[is.nan(D.estimate.table$Estimate)] <- 0
    D.estimate.table$se[is.nan(D.estimate.table$se)] <- 0
    D.estimate.table$cv[is.nan(D.estimate.table$cv)] <- 0
    D.estimate.table$lcl[is.nan(D.estimate.table$lcl)] <- 0
    D.estimate.table$ucl[is.nan(D.estimate.table$ucl)] <- 0
    
    # Return list depending on value of se
    # change to set missing values to 0
    # jll 6/30/06; dropped restriction that numregions > 1 on sending vc back
    if(se){
      cormat <- result$vc/(result$estimate.table$se %o% result$estimate.table$se)
      cormat[is.nan(cormat)] <- 0
      result <- list(bysample=bysample.table, summary = summary.table,
                     N=result$estimate.table, D=D.estimate.table, 
                     average.p=average.p, cormat = cormat,
                     vc=list(total=result$vc,detection=result$vc1,er=result$vc2),
                     Nhat.by.sample=Nhat.by.sample)
    }else{
      result <- list(bysample=bysample.table,summary = summary.table,
                     N = estimate.table,D = D.estimate.table, average.p=average.p,
                     Nhat.by.sample=Nhat.by.sample)
    }
    return(result)
  }
  
  ###Start of dht function
  
  # Code additions by jll 18-Nov-04; the following allows for a subset
  # statement to be added to create obs.table from model data rather than
  # creating obs.table separately. This only works if the data contain the
  # Sample.Label and Region.Label fields.
  point <- model$meta.data$point
  objects <- as.numeric(names(model$fitted))
  if(is.null(obs.table)){
    data <- model$data
    if("observer"%in%names(data)){
      # jll 3 Sept 2014 if dual observer I added code to use observer 1 only
      # or it was doubling sample size
      data <- data[data$observer==1,]
    }
    if("Sample.Label" %in% names(data) & "Region.Label" %in% names(data)){
      if(is.null(substitute(subset))){
        obs.table <- data[,c("object","Sample.Label","Region.Label")]
      }else{
        select <- data[eval(substitute(subset),envir=data),]
        obs.table <- select[,c("object","Sample.Label","Region.Label")]
      }
      obs.table <- obs.table[obs.table$object %in% objects,]
    }else{
      stop("Must specify obs.table because Sample.Label and/or Region.Label fields not contained in data")
    }
  }
  
  # Extract relevant fields from Region and Sample tables; jll 4 May 07;
  region.table <- region.table[,c("Region.Label","Area")]
  sample.table <- sample.table[,c("Region.Label","Sample.Label","Effort")]
  
  # Make sure input data labels are factors
  region.table$Region.Label <- factor(region.table$Region.Label)
  sample.table$Region.Label <- factor(sample.table$Region.Label,
                                      levels=levels(region.table$Region.Label))
  obs.table$Region.Label <- factor(obs.table$Region.Label,
                                   levels=levels(region.table$Region.Label))
  sample.table$Sample.Label <- factor(sample.table$Sample.Label)
  obs.table$Sample.Label <- factor(obs.table$Sample.Label,
                                   levels=levels(sample.table$Sample.Label))
  
  # Assign default values to options
  options <- assign.default.values(options, pdelta = 0.001, varflag = 2,
                                   convert.units = 1, ervar="R2")
  
  # Convert width value
  width <- model$meta.data$width * options$convert.units
  
  # If area is zero for all regions reset to the area of the covered region
  DensityOnly <- FALSE
  if(sum(region.table$Area)==0){
    DensityOnly <- TRUE
    # cat("Warning: Area for regions is zero. They have been set to area of covered region(strips), \nso N is for covered region.",
    #     "However, standard errors will not match \nprevious covered region SE because it includes spatial variation\n")
    Effort.by.region <- by(sample.table$Effort, sample.table$Region.Label,sum)
    region.table$Area <- ifelse(point,
                                pi*as.vector(Effort.by.region)*width^2,
                                2*as.vector(Effort.by.region)*width)
  }
  
  # Create obs/samples structures
  vs <- create.varstructure(model, region.table, sample.table, obs.table)
  samples <- vs$samples
  obs <- vs$obs
  region.table <- vs$region
  
  # handle subset feature when labels are also in data
  if(!is.null(obs$Region.Label.x)){
    obs$Region.Label <- obs$Region.Label.x
    obs$Sample.Label <- obs$Sample.Label.x
    obs$Region.Label.x <- NULL
    obs$Sample.Label.x <- NULL
    obs$Region.Label.y <- NULL
    obs$Sample.Label.y <- NULL
  }
  
  # Merge with fitted values
  pdot <- model$fitted
  obs <- merge(obs,data.frame(object=objects,pdot=pdot))
  
  # If clustered population create tables for clusters and individuals and
  # an expected S table otherwise just tables for individuals in an
  # unclustered popn
  if(!is.null(obs$size)){
    clusters <- tables.dht(TRUE)
    individuals <- tables.dht(FALSE)
    Expected.S <- individuals$N$Estimate/clusters$N$Estimate
    
    # This computes the se(E(s)).  It essentially uses 3.37 from Ads but in
    # place of using 3.25, 3.34 and 3.38, it uses 3.27,3.35 and an equivalent
    # cov replacement term for 3.38.  This uses line to line variability
    # whereas the other formula measure the variance of E(s) within the lines 
    # and it goes to zero as p approaches 1.
    if(se & options$varflag!=1){
      if(options$varflag==2){
        numRegions <- length(unique(samples$Region.Label))
        cov.Nc.Ncs <- rep(0,numRegions)
        scale <- clusters$summary$Area/clusters$summary$CoveredArea
        
        for(i in 1:numRegions){
          c.stratum.data <- clusters$Nhat.by.sample[
            as.character(clusters$Nhat.by.sample$Region.Label) ==
              as.character(region.table$Region.Label[i]), ]
          
          i.stratum.data <- individuals$Nhat.by.sample[
            as.character(individuals$Nhat.by.sample$Region.Label) ==
              as.character(region.table$Region.Label[i]), ]
          
          Li <- sum(c.stratum.data$Effort.x)
          cov.Nc.Ncs[i] <- covn(c.stratum.data$Effort.x/(scale[i]*Li),
                                c.stratum.data$Nhat/scale[i],
                                i.stratum.data$Nhat/scale[i],
                                options$ervar)
        }
      }else{
        cov.Nc.Ncs <- as.vector(by(obs$size*(1 - obs$pdot)/obs$pdot^2,
                                   obs$Region.Label, sum))
        cov.Nc.Ncs[is.na(cov.Nc.Ncs)] <- 0
      }
      
      cov.Nc.Ncs[is.nan(cov.Nc.Ncs)] <- 0
      cov.Nc.Ncs <- c(cov.Nc.Ncs,sum(cov.Nc.Ncs))
      cov.Nc.Ncs <- cov.Nc.Ncs+diag(t(clusters$vc$detection$partial)%*%
                                      solvecov(model$hessian)$inv%*%
                                      individuals$vc$detection$partial)
      se.Expected.S <- clusters$N$cv^2 + individuals$N$cv^2 -
        2*cov.Nc.Ncs/(individuals$N$Estimate*clusters$N$Estimate)
      Expected.S[is.nan(Expected.S)] <- 0
      se.Expected.S[se.Expected.S<=0 | is.nan(se.Expected.S)] <- 0
      se.Expected.S <- Expected.S*sqrt(se.Expected.S)
      Expected.S <- data.frame(Region=clusters$N$Label,
                               Expected.S=as.vector(Expected.S),
                               se.Expected.S=as.vector(se.Expected.S))
    }else{
      Expected.S <- data.frame(Region=clusters$N$Label,
                               Expected.S=as.vector(Expected.S))
    }
    
    if(DensityOnly){
      clusters$N <- NULL
      individuals$N <- NULL
    }
    
    result <- list(clusters=clusters,
                   individuals=individuals,
                   Expected.S=as.vector(Expected.S))
  }else{
    individuals <- tables.dht(TRUE)
    if(DensityOnly){
      individuals$N <- NULL
    }
    result <- list(individuals=individuals)
  }
  
  class(result) <- "dht"
  return(result)
}

survey.region.dht <- function(Nhat.by.sample, samples, width, point){
  #  Compute effort in each region and the area in the covered region
  Effort.by.region <- by(samples$Effort,samples$Region.Label,sum)
  if(point){
    CoveredArea <- pi*as.vector(Effort.by.region)*width^2
  }else{
    CoveredArea <- as.vector(Effort.by.region)*width ####removed *2 to make one sided####
  }
  
  # Scale up abundance in covered region to the survey region
  # unless no areas given
  Nhat.by.sample <- merge(Nhat.by.sample,
                          data.frame(Region.Label=names(Effort.by.region),
                                     CoveredArea=CoveredArea,
                                     Effort=as.vector(Effort.by.region)),
                          by.x="Region.Label",
                          by.y="Region.Label",
                          all.x=TRUE)
  Nhat.by.sample$Nhat <- Nhat.by.sample$Nhat*Nhat.by.sample$Area/
    Nhat.by.sample$CoveredArea
  return(Nhat.by.sample)
}

#' Variance and confidence intervals for density and abundance estimates
#'
#' Computes standard error, cv, and log-normal confidence intervals for
#' abundance and density within each region (if any) and for the total of all
#' the regions.  It also produces the correlation matrix for regional and total
#' estimates.
#'
#' The variance has two components: 1) variation due to uncertanity from
#' estimation of the detection function and 2) variation in abundance due to
#' random sample selection.  The first component is computed using a delta
#' method estimate of variance (\code{\link{DeltaMethod}} (Huggins 1989, 1991,
#' Borchers et al. 1998) in which the first derivatives of the abundance
#' estimator with respect to the parameters in the detection function are
#' computed numerically.  The second component can be computed in one of three
#' ways as set by the option \code{varflag} with values 0,1,2.
#'
#' A value of 0 is to use a binomial variance for the number of observations
#' and it is only useful if the sampled region is the survey region and the
#' objects are not clustered which will not occur very often. If covered region
#' is less than the survey region the variance estimator is scaled up but it
#' will be a poor estimator and the confidence interval will likely not achieve
#' the nominal level.
#'
#' A value of 1 uses the variance for the encounter rate of (Fewster et al.
#' (2009), estimator R2 (which has been shown to have better properties than
#' the previous default of Buckland et al. 2001 pg 78-79)).  If
#' \code{group=FALSE} the variance of the mean group size is also included.
#' This variance estimator is not appropriate if \code{size} or a derivative of
#' \code{size} is used in the any of the detection function models.
#'
#' In general if any covariates are used in the models, the default option 2 is
#' preferable.  It uses a variance estimator based on that suggested by Innes
#' et al. (2002) which used the formula for the variance ecounter rate but
#' replaces the number of observations per sample with the estimated abundance
#' per sample.  The difference between the version used here and that in Innes
#' et al. (2002) is that Innes et al. use an estimator with form similar to
#' that of Buckland et al. (2001), while the estimator here uses a form based
#' on Fewster et al. (2009, estimator R2).
#'
#' For more on encounter rate variance estimation, see \link{varn}.
#'
#' Exceptions to the above occur if there is only one sample in a stratum. In
#' that case it uses Poisson assumption (var(x)=x) and it assumes a known
#' variance so z=1.96 is used for critical value. In all other cases the
#' degrees of freedom for the t-distribution assumed for the log(abundance) or
#' log(density) is based on the Satterthwaite approximation (Buckland et al.
#' 2001 pg 90) for the degrees of freedom (df).  The df are weighted by the
#' squared cv in combining the two sources of variation because of the assumed
#' log-normal distribution because the components are multiplicative.  For
#' combining df for the sampling variance across regions they are weighted by
#' the variance because it is a sum across regions.
#'
#' A non-zero correlation between regional estimates can occur from using a
#' common detection function across regions.  This is reflected in the
#' correlation matrix of the regional and total estimates which is given in the
#' value list.  It is only needed if subtotals of regional estimates are
#' needed.
#'
#' @param model ddf model object
#' @param region.table table of region values
#' @param samples table of samples(replicates)
#' @param obs table of observations
#' @param options list of options that can be set (see \code{\link{dht}})
#' @param numRegions number of regions
#' @param estimate.table table of estimate values
#' @param Nhat.by.sample estimated abundances by sample
#' @export
#' @return List with 2 elements: \item{estimate.table}{completed table with se,
#'   cv and confidence limits} \item{vc }{correlation matrix of estimates}
#' @note This function is called by \code{dht} and it is not expected that the
#'   user will call this function directly but it is documented here for
#'   completeness and for anyone expanding the code or using this function in
#'   their own code
#' @author Jeff Laake
#' @seealso \code{\link{dht}}, \code{\link{print.dht}}
#' @references see \code{\link{dht}}
#' @keywords utility
#' @importFrom stats qnorm qt var
dht.se <- function(model, region.table, samples, obs, options, numRegions,
                   estimate.table, Nhat.by.sample){
  #  Functions Used:  DeltaMethod, dht.deriv (in DeltaMethod), varn
  
  # Define function: compute.df
  compute.df<- function(k,type){
    if(type=="O1" | type=="O2"| type=="O3"){
      H.O <- k - 1
      k.h.O <- rep(2, H.O)
      df <- sum(k.h.O - 1)
    }else{
      if(type=="S1" | type=="S2"){
        H.S <- floor(k/2)
        k.h.S <- rep(2, H.S)
        if(k %% 2 > 0) k.h.S[H.S] <- 3
        df <- sum(k.h.S - 1)
      }else{
        df <- k-1
      }
    }
    return(df)
  }
  
  # First compute variance component due to estimation of detection function
  # parameters. This uses the delta method and produces a v-c matrix if more
  # than one strata
  if(!is.null(model$par)){
    vcov <- solvecov(model$hessian)$inv
    vc1.list <- DeltaMethod(model$par, dht.deriv, vcov, options$pdelta,
                            model = model, samples = samples, obs = obs, options = options)
    vc1 <- vc1.list$variance
  }else{
    vc1.list <- list(variance=0)
    vc1 <- 0
  }
  
  # Next compute the component due to sampling of both lines and of the
  # detection process itself
  # There are 3 different options here:
  #  1) varflag=0; Binomial variance of detection process - only applicable if
  #   survey region=covered region although it will scale up but it would be
  #   a poor estimator
  #  2) varflag=1; delta method, with varn based on Fewster et al (2009)
  #   estimator R2 (var(n/L))
  #  3) varflag=2; Innes et al variance estimator (var(N/L), except changed to
  #   resemble the form of estimator R2 of Fewster et al (2009))
  # Exceptions to the above occur if there is only one sample in a stratum.
  #   In that case it uses Poisson approximation.
  
  scale <- region.table$Area/region.table$CoveredArea
  
  # If no areas were given or varflag=0 use approach #1 (varflag=0)
  # Note: vc2 is of proper dimension because Region.Label for obs is setup
  # with all levels of the Region.Label from the region.table.
  if(sum(region.table$Area) == 0 | options$varflag == 0){
    if(options$group){
      vc2 <- by((1 - obs$pdot)/obs$pdot^2, obs$Region.Label,sum)
    }else{
      vc2 <- by(obs$size^2 * (1 - obs$pdot)/obs$pdot^2, obs$Region.Label, sum)
    }
    # set missing value to 0
    vc2[is.na(vc2)] <- 0
    
    if(sum(region.table$Area) != 0){
      vc2 <- vc2 * scale[1:numRegions]^2
    }
  }else{
    # Otherwise compute variance for varflag=1 or 2
    vc2 <- rep(0, numRegions)
    # 26 jan 06 jll; changed to use object rather than distance; also
    # overwrites existing n because that can be sum(size) rather than count
    nobs <- tapply(obs$object, obs$Label, length)
    nobs <- data.frame(Label = names(nobs),
                       n = as.vector(nobs)[!is.na(nobs)])
    Nhat.by.sample$n <- NULL
    # when there are no sighings
    if(nrow(nobs) > 0){
      Nhat.by.sample <- merge(Nhat.by.sample, nobs, by.x = "Label",
                              by.y = "Label", all.x = TRUE)
      Nhat.by.sample$n[is.na(Nhat.by.sample$n)] <- 0
    }else{
      Nhat.by.sample <- cbind(Nhat.by.sample, n = rep(0,nrow(Nhat.by.sample)))
    }
    
    # Compute number of lines per region for df calculation
    if(numRegions > 1){
      estimate.table$k <- c(as.vector(tapply(samples$Effort,
                                             samples$Region.Label, length)), 0)
      estimate.table$k[numRegions + 1] <- sum(estimate.table$k)
    }else{
      estimate.table$k <- as.vector(tapply(samples$Effort,
                                           samples$Region.Label, length))
    }
    
    # If individual abundance being computed, calculate mean and variance
    # of mean for group size.
    if(!options$group){
      if(length(obs$size) > 0){
        vars <- by(obs$size, obs$Region.Label, var)/
          by(obs$size, obs$Region.Label, length)
        sbar <- by(obs$size, obs$Region.Label, mean)
        sobs <- data.frame(Region.Label = names(sbar),
                           vars         = as.vector(vars),
                           sbar         = as.vector(sbar))
      }else{
        sobs = data.frame(Region.Label=levels(obs$Region.Label),
                          vars=rep(NA,length(levels(obs$Region.Label))),
                          sbar=rep(NA,length(levels(obs$Region.Label))))
      }
      Nhat.by.sample <- merge(Nhat.by.sample, sobs, by.x = "Region.Label",
                              by.y = "Region.Label", all.x = TRUE)
      Nhat.by.sample$sbar[is.na(Nhat.by.sample$sbar)] <- 0
      Nhat.by.sample$vars[is.na(Nhat.by.sample$vars)] <- 0
    }else{
      # If group abundance is being estimated, set mean=1, var=0
      Nhat.by.sample$sbar <- rep(1, dim(Nhat.by.sample)[1])
      Nhat.by.sample$vars <- rep(0, dim(Nhat.by.sample)[1])
    }
    
    # sort Nhat.by.sample by Region.Label and Sample.Label
    Nhat.by.sample <- Nhat.by.sample[order(Nhat.by.sample$Region.Label,Nhat.by.sample$Sample.Label),]
    
    # Loop over each region and compute each variance;
    # jll 11/11/04 - changes made in following code using
    # Effort.x (effort per line) rather than previous errant code
    # that used Effort.y (effort per region)
    for(i in 1:numRegions){
      stratum.data <- Nhat.by.sample[as.character(Nhat.by.sample$Region.Label)==
                                       as.character(region.table$Region[i]), ]
      Ni <- sum(stratum.data$Nhat)
      Li <- sum(stratum.data$Effort.x)
      sbar <- stratum.data$sbar[1]
      vars <- stratum.data$vars[1]
      
      if (options$group) vars <- 0
      
      if(length(stratum.data$Effort.y) == 1){
        if (options$varflag == 1){
          vc2[i] <- Ni^2 * (1/stratum.data$n + vars/sbar^2)
        }else{
          vc2[i] <- Ni^2 * (1/Ni + vars/sbar^2)
        }
      }else if (options$varflag == 1){
        vc2[i] <- (Ni * Li)^2 * varn(stratum.data$Effort.x,
                                     stratum.data$n,type=options$ervar)/
          sum(stratum.data$n)^2 + Ni^2 * vars/sbar^2
      }else{
        if(options$varflag==2){
          vc2[i] <- varn(stratum.data$Effort.x/(scale[i] * Li),
                         stratum.data$Nhat/scale[i],type=options$ervar)
        }else{
          vc2[i] <- varn(stratum.data$Effort.x/(scale[i] * Li),
                         stratum.data$Nhat/scale[i], type=options$ervar)
        }
      }
    }
  }
  
  vc2[is.nan(vc2)] <- 0
  
  # Construct v-c matrix for encounter rate variance given computed ps
  # The cov between regional estimate and total estimate is simply var for
  #  regional estimate.
  # Assumes no cov between regions due to independent sample selection.
  if(numRegions > 1){
    v2 <- vc2
    vc2 <- diag(c(vc2, sum(vc2)))
    vc2[1:numRegions, (numRegions + 1)] <- v2
    vc2[(numRegions + 1), 1:numRegions] <- v2
  }else if (length(vc2) > 1){
    vc2 <- diag(vc2)
  }else{
    vc2 <- as.matrix(vc2)
  }
  
  vc <- vc1 + vc2
  
  # deal with missing values and 0 estimates.
  estimate.table$se <- sqrt(diag(vc))
  estimate.table$se[is.nan(estimate.table$se)] <- 0
  estimate.table$cv <- estimate.table$se/estimate.table$Estimate
  estimate.table$cv[is.nan(estimate.table$cv)] <- 0
  
  # work out the confidence intervals
  # if the options$ci.width is set, then use that, else default to
  # 95% CI
  if(is.null(options$ci.width)){
    ci.width <- 0.025
  }else{
    ci.width <- (1-options$ci.width)/2
  }
  
  # Use satterthwaite approx for df and log-normal distribution for
  # 95% intervals
  if(options$varflag != 0){
    # set df from replicate lines to a min of 1 which avoids divide by zero
    # Following 2 lines added and references to estimate.table$k changed to df
    df <- estimate.table$k
    df <- sapply(df,compute.df,type=options$ervar)
    df[df<1] <- 1
    
    if(is.na(vc1) || all(vc1==0)){
      estimate.table$df <- df
    }else{
      estimate.table$df <- estimate.table$cv^4/((diag(vc1)/
                                                   estimate.table$Estimate^2)^2/(length(model$fitted) -
                                                                                   length(model$par)) +
                                                  (diag(vc2)/estimate.table$Estimate^2)^2/df)
    }
    
    # compute proper satterthwaite
    # df for total estimate assuming sum of indep region estimates; uses
    # variances instead of cv's because it is a sum of means for encounter
    # rate portion of variance (df.total)
    if(numRegions>1){
      df.total <- (diag(vc2)[numRegions+1])^2/
        sum((diag(vc2)^2/df)[1:numRegions])
      if(all(vc1==0)){
        estimate.table$df[numRegions+1] <- df.total
      }else{
        estimate.table$df[numRegions+1] <- estimate.table$cv[numRegions+1]^4 /
          ((diag(vc1)[numRegions+1]/
              estimate.table$Estimate[numRegions+1]^2)^2/
             (length(model$fitted)-length(model$par))
           + (diag(vc2)[numRegions+1]/
                estimate.table$Estimate[numRegions+1]^2)^2/df.total)
      }
    }
    
    estimate.table$df[estimate.table$df < 1 &estimate.table$df >0] <- 1
    cvalue <- exp((abs(qt(ci.width, estimate.table$df)) *
                     sqrt(log(1 + estimate.table$cv^2))))
  }else{
    # intervals for varflag=0; sets df=0
    # and uses normal approximation
    cvalue <- exp((abs(qnorm(ci.width)) * sqrt(log(1 + estimate.table$cv^2))))
    estimate.table$df <- rep(0,dim(estimate.table)[1])
  }
  
  # deal with missing values and divide by 0 issues
  estimate.table$df[is.nan(estimate.table$df)] <- 0
  estimate.table$lcl  <-  estimate.table$Estimate/cvalue
  estimate.table$lcl[is.nan(estimate.table$lcl)] <- 0
  estimate.table$ucl  <-  estimate.table$Estimate * cvalue
  estimate.table$ucl[is.nan(estimate.table$ucl)] <- 0
  estimate.table$k  <-  NULL
  
  print(c(vc, vc1, vc2))
  
  return(list(estimate.table = estimate.table,
              vc             = vc,
              vc1            = vc1.list,
              vc2            = vc2 ))
}

#' Assign default values to list elements that have not been already assigned
#'
#' Assigns default values for \code{argument} in list \code{x} from
#' \code{argument=value} pairs in \dots{} if \code{x$argument} doesn't already
#' exist
#'
#' @param x generic list
#' @param \dots unspecified list of argument=value pairs that are used to
#'   assign values
#' @return x - list with filled values
#' @author Jeff Laake
#' @keywords ~utility
assign.default.values <- function(x,...){
  args <- list(...)
  notset <- !names(args)%in%names(x)
  x <- c(x,args[notset])
  return(x)
}

#' Creates structures needed to compute abundance and variance
#'
#' Creates samples and obs dataframes used to compute abundance and its
#' variance based on a structure of geographic regions and samples within each
#' region.  The intent is to generalize this routine to work with other
#' sampling structures.
#'
#' The function performs the following tasks: 1)tests to make sure that region
#' labels are unique, 2) merges sample and region tables into a samples table
#' and issue a warning if not all samples were used, 3) if some regions have no
#' samples or if some values of Area were not valid areas given then issue
#' error and stop, then an error is given and the code stops, 4) creates a
#' unique region/sample label in samples and in obs, 5) merges observations
#' with sample and issues a warning if not all observations were used, 6) sorts
#' regions by its label and merges the values with the predictions from the
#' fitted model based on the object number and limits it to the data that is
#' appropriate for the fitted detection function.
#'
#' @param model fitted ddf object
#' @param region region table
#' @param sample sample table
#' @param obs table of object #'s and links to sample and region table
#' @return List with 2 elements: \item{samples }{merged dataframe containing
#'   region and sample info - one record per sample} \item{obs}{merged
#'   observation data and links to region and samples}
#' @note Internal function called by \code{\link{dht}}
#' @author Jeff Laake
#' @keywords utility
create.varstructure <- function(model,region,sample,obs){
  
  # Test to make sure that region labels are unique
  if(length(unique(region$Region.Label))!=length(region$Region.Label)){
    stop("Region labels must be unique")
  }
  
  # Merge sample and region tables into samples; warn if not all samples used
  samples <- merge(region,sample,by.x="Region.Label",all.x=TRUE,all.y=TRUE)
  if(any(is.na(samples$Area))){
    warning("Some samples not included in the analysis")
  }
  
  # Test to make sure that sample labels are unique within region
  if(dim(samples)[1] != dim(unique(data.frame(region=samples$Region.Label,
                                              sample=samples$Sample.Label)))[1]){
    stop("Sample labels must be unique within region")
  }
  
  # Merge again but don't use all.y which ignores samples not used
  samples <- merge(region,sample,by.x="Region.Label",all.x=TRUE)
  samples <- samples[order(samples$Region.Label,samples$Sample.Label),]
  
  # If some regions have no samples then issue error and stop; also
  # if invalid areas given then issue error and stop
  if(any(is.na(samples$Sample.Label))){
    stop(paste("Following regions have no samples - ",
               paste(samples$Region.Label[is.na(samples$Sample.Label)],collapse=",")))
  }
  
  if(any(is.na(region$Area)) | any(!is.numeric(region$Area))){
    stop("Invalid or missing Area values for regions\n")
  }
  
  # Create a unique region/sample label in samples and in obs
  samples$Label <- paste(samples$Region.Label,samples$Sample.Label,sep="")
  obs$Label <- paste(obs$Region.Label,obs$Sample.Label,sep="")
  
  # we only want the following columns from obs:
  #  Label, object, Region.Label, Sample.Label
  # so remove everything else so we don't end up with .x and .ys in the 
  # merges that follow...
  obs <- obs[,c("Label", "object", "Region.Label", "Sample.Label")]
  
  # Merge observations with sample; warn if not all observations used
  data <- merge(obs,samples,by.x="Label",by.y="Label",all.x=TRUE,sort=FALSE)
  
  if(any(is.na(data$Region.Label.y))){
    warning("Some observations not included in the analysis")
  }
  
  data <- data[!is.na(data$Region.Label.y),]
  data$Region.Label.y <- NULL
  data$Sample.Label.y <- NULL
  data$Sample.Label <- data$Sample.Label.x
  data$Region.Label <- data$Region.Label.x
  data$Region.Label.x <- NULL
  data$Sample.Label.x <- NULL
  
  # Sort regions by label
  region <- region[order(region$Region.Label),]
  
  # Merge with data from model and limit to data appropriate for method
  data <- merge(data,model$data,by.x="object",by.y="object",sort=FALSE)
  
  # observer =1 to avoid problems with merge; this forces abundance 
  #  to always be estimated using observer 1 as the primary
  obs <- 1
  if(!model$method %in% c("io","io.fi","rem","rem.fi")){
    if(model$method!="ds"){
      data <- data[data$observer==obs,]
      data <- data[data$detected==1,]
    }
  }else{
    data <- data[data$observer==obs,]
  }
  
  # Return vectors and lists for computation
  return(list(samples=samples,
              obs=data,
              region=region))
}

#' Covered region estimate of abundance from Horvitz-Thompson-like estimator
#'
#' Computes H-T abundance within covered region by sample.
#'
#'
#' @param obs observations table
#' @param samples samples table
#' @param group if TRUE compute abundance of group otherwise abundance of
#'   individuals
#' @return Nhat.by.sample - dataframe of abundance by sample
#' @note Internal function called by \code{\link{dht}} and related functions
#' @author Jeff Laake
#' @keywords utility
covered.region.dht <- function(obs, samples, group){
  
  # Compute abundance in covered region depending on value of
  # group = TRUE (do group abundance); F(do individual abundance)
  # if there are observations of this species
  if(nrow(obs) > 0){
    Nhats <- compute.Nht(obs$pdot, group, obs$size)
    
    # Sum abundances by sample within region
    Nhats <- by(Nhats, obs$Label, sum)
    
    # Sum observations by sample within region
    if(group){
      sum.obs <- by(obs$object, obs$Label, length)
    
    }else{
      sum.obs <- by(obs$size, obs$Label, sum)
    }
    
    # Merge with samples
    num.obs <- data.frame(Label = names(sum.obs), n = as.vector(sum.obs))
    Nhats <- data.frame(Label = names(Nhats), Nhat = as.vector(Nhats))
    Nhat.by.sample <- merge(samples, Nhats,by.x = "Label", by.y = "Label", all.x = TRUE)
    Nhat.by.sample <- merge(Nhat.by.sample, num.obs, by.x = "Label", by.y = "Label", all.x = TRUE)
    Nhat.by.sample$Nhat[is.na(Nhat.by.sample$Nhat)] <- 0
    Nhat.by.sample$n[is.na(Nhat.by.sample$n)] <- 0
    
    # Create Nhat.by.sample in the case where there were no sightings
  }else{
    Nhat.by.sample <- cbind(samples,
                            Nhat = rep(0,nrow(samples)),
                            n = rep(0,nrow(samples)))
  }
  return(Nhat.by.sample)
}

compute.Nht <- function(pdot,group=TRUE,size=NULL){
  if(group){
    return(1/pdot)
  }else{
    if(!is.null(size)){
      if(!any(is.na(size))){
        return(size/pdot)
      }else{
        stop("One or more missing group sizes in observations")
      }
    }else{
      stop("Missing group sizes in observations")
    }
  }
}

#' Compute empirical variance of encounter rate
#'
#' Computes one of a series of possible variance estimates for the observed
#' encounter rate for a set of sample measurements (e.g., line lengths) and
#' number of observations per sample.
#'
#' The choice of type follows the notation of Fewster et al. (2009) in that
#' there are 8 choices of encounter rate variance that can be computed:
#'
#' \describe{
#' \item{\code{R2}}{random line placement with unequal line lengths (design-assisted estimator)}
#' \item{\code{R3}}{random line placement, model-assisted estimator, based on true contagion process}
#' \item{\code{R4}}{random line placement, model-assisted estimator, based on apparent contagion process}
#' \item{\code{S1}}{systematic line placement, post-stratification with no strata overlap}
#' \item{\code{S2}}{systematic line placement, post-stratification with no strata overlap, variances weighted by line length per stratum}
#' \item{\code{O1}}{systematic line placement, post-stratification with overlapping strata (akin to S1)}
#' \item{\code{O2}}{systematic line placement, post-stratification with overlapping strata (weighted by line length per stratum, akin to S2)}
#' \item{\code{O3}}{systematic line placement, post-stratification with overlapping strata, model-assisted estimator with trend in encounter rate with line length}}
#'
#' Default value is R2, shown in Fewster et al. (2009) to have good performance
#' for completely random designs.  For systematic parallel line transect
#' designs, Fewster et al. recommend O2.
#'
#' For the systematic estimators, pairs are assigned in the order they are
#' given in the \code{lengths} and \code{groups} vectors.
#'
#' @usage   varn(lvec,nvec,type)
#'
#'          covn(lvec, groups1, groups2, type)
#' @aliases varn covn
#' @param lvec vector of sample measurements (e.g., line lengths)
#' @param nvec vector of number observed
#' @param groups1 vector of number of groups observed
#' @param groups2 vector of number of individuals observed
#' @param type choice of variance estimator to use for encounter rate
#' @return Variance of encounter rate as defined by arguments
#' @note This function is also used with different calling arguments to compute
#'   Innes et al variance of the estimated abundances/length rather than
#'   observation encounter rate. The function covn is probably only valid for
#'   R3 and R2.  Currently, the R2 form is used for all types other than R3.
#' @author Jeff Laake
#' @references Fewster, R.M., S.T. Buckland, K.P. Burnham, D.L. Borchers, P.E.
#'   Jupp, J.L. Laake and L. Thomas. 2009. Estimating the encounter rate
#'   variance in distance sampling. Biometrics 65: 225-236.
#' @keywords utility
varn <- function(lvec,nvec,type){
  #  Function courtesy of Rachel Fewster with a few minor changes
  
  ntot <- sum(nvec)
  L <- sum(lvec)
  k <- length(lvec)
  
  ## Go through the estimators one by one.
  ## R2, R3, R4, S1, S2, O1, O2, O3
  if(!(type %in% c("R2","R3","R4","S1","S2","O1","O2","O3")))
    stop (paste("Encounter rate variance type '",type,"' is not recognized.",sep=""))
  
  ## First the estimators based on the assumption of a random sample
  ## of lines: R2, R3, R4:
  
  ## Estimator R2: var.R2, sd.R2
  if(type=="R2"){
    var.R2 <- (k * sum(lvec^2 * (nvec/lvec - ntot/L)^2))/(L^2 * (k -1))
    return(var.R2)
  }
  
  ## Estimator R3: var.R3, sd.R3
  if(type=="R3"){
    var.R3 <- 1/(L * (k - 1)) * sum(lvec * (nvec/lvec - ntot/L)^2)
    return(var.R3)
  }
  
  ## Estimator R4: var.R4, sd.R4
  ## Using approximation to Negbin variance (the apparent contagion model)
  ## using phi = 2 + eps
  if(type=="R4"){
    if(all(lvec==mean(lvec))){
      phi <- (2-2/k)/(1-2/k)
    }else{
      S <- sum(lvec^2)
      C <- sum(lvec^3)
      logvec <- log(lvec)
      D1 <- sum(lvec * logvec)
      D2 <- sum(lvec^2 * logvec)
      D3 <- sum(lvec^3 * logvec)
      eps.top <- 2 * (S^2 - L * C)
      eps.bottom <- L * S * D1 - 2 * S * D2 + 2 * L * D3 - L^2 * D2
      eps <- eps.top/eps.bottom
      phi <- 2 + eps
    }
    alpha <- 1/sum(lvec^phi * (L/lvec - 1))
    var.R4 <- alpha * sum(lvec^phi * (nvec/lvec - ntot/L)^2)
    return(var.R4)
  }
  
  ## Now the stratified estimators with non-overlapping strata:
  ## S1 and S2:
  
  ## First group the lines into strata, so that all strata have
  ## two lines but if the last stratum has three if necessary:
  H <- floor(k/2)
  k.h <- rep(2, H)
  if(k %% 2 > 0){
    k.h[H] <- 3
  }
  end.strat <- cumsum(k.h)
  begin.strat <- cumsum(k.h) - k.h + 1
  
  ## Estimators S1 and S2: var.S1, sd.S1 ; var.S2, sd.S2
  if(type=="S1" | type=="S2"){
    sum.S1 <- 0
    sum.S2 <- 0
    for(h in 1:H){
      nvec.strat <- nvec[begin.strat[h]:end.strat[h]]
      lvec.strat <- lvec[begin.strat[h]:end.strat[h]]
      nbar.strat <- mean(nvec.strat)
      lbar.strat <- mean(lvec.strat)
      
      ## S1 calculations:
      inner.strat.S1 <- sum((nvec.strat - nbar.strat - (ntot/L) *
                               (lvec.strat - lbar.strat))^2)
      sum.S1 <- sum.S1 + k.h[h]/(k.h[h] - 1) * inner.strat.S1
      
      ## S2 calculations: note that we use estimator R2 within
      ## each stratum:
      L.strat <- sum(lvec.strat)
      var.strat.S2 <- k.h[h]/(L.strat^2 * (k.h[h] - 1)) * sum(
        lvec.strat^2 * (nvec.strat/lvec.strat - nbar.strat/lbar.strat)^2)
      sum.S2 <- sum.S2 + L.strat^2 * var.strat.S2
    }
    if(type=="S1"){
      var.S1 <- sum.S1/L^2
      return(var.S1)
    }else{
      var.S2 <- sum.S2/L^2
      return(var.S2)
    }
  }
  
  ## Now the stratified estimators with overlapping strata:
  ## O1, O2, O3:
  lvec.1 <- lvec[ - k]
  lvec.2 <- lvec[-1]
  nvec.1 <- nvec[ - k]
  nvec.2 <- nvec[-1]
  ervec.1 <- nvec.1/lvec.1
  ervec.2 <- nvec.2/lvec.2
  
  ## Estimator O1: var.O1, sd.O1
  if(type=="O1"){
    overlap.varterm <- (nvec.1 - nvec.2 - ntot/L * (lvec.1 - lvec.2))^2
    var.O1 <- k/(2 * L^2 * (k - 1)) * sum(overlap.varterm)
    return(var.O1)
  }
  
  ## Estimator O2: var.O2, sd.O2
  if(type=="O2"){
    V.overlap.R2 <- ((lvec.1 * lvec.2)/(lvec.1 + lvec.2))^2 * (ervec.1-ervec.2)^2
    var.O2 <- (2 * k)/(L^2 * (k - 1)) * sum(V.overlap.R2)
    return(var.O2)
  }
  
  ## Estimator O3: var.O3, sd.O3
  if(type=="O3"){
    V.overlap.R3 <- ((lvec.1 * lvec.2)/(lvec.1 + lvec.2)) * (ervec.1 -ervec.2)^2
    var.O3 <- 1/(L * (k - 1)) * sum(V.overlap.R3)
    return(var.O3)
  }
}

covn <- function (lvec, groups1, groups2, type){
  # covn - computes covariance of encounter rate of clusters and individuals as
  #        called from dht.  modeled after varn.
  L <- sum(lvec)
  n1 <- sum(groups1)
  er1 <- n1/L
  n2 <- sum(groups2)
  er2 <- n2/L
  
  k=length(lvec)
  if(type=="R3"){
    varer <- sum(lvec * (groups1/lvec - er1)*(groups2/lvec-er2))/(L*(k-1))
  }else{
    varer <- k*sum(lvec^2 * (groups1/lvec - er1)*(groups2/lvec-er2))/(L^2*(k-1))
  }
  
  return(varer)
}

# solvecov code was taken from package fpc: Christian
# Hennig chrish@@stats.ucl.ac.uk http://www.homepages.ucl.ac.uk/~ucakche/
solvecov <- function (m, cmax = 1e+10){
  options(show.error.messages = FALSE)
  covinv <- try(solve(m))
  if(class(covinv) != "try-error"){
    coll = FALSE
  }else{
    p <- nrow(m)
    cove <- eigen(m, symmetric = TRUE)
    coll <- TRUE
    if(min(cove$values) < 1/cmax) {
      covewi <- diag(p)
      for (i in 1:p){
        if (cove$values[i] < 1/cmax){
          covewi[i, i] <- cmax
        }else{
          covewi[i, i] <- 1/cove$values[i]
        }
      }
    }else{
      covewi <- diag(1/cove$values, nrow = length(cove$values))
    }
    covinv <- cove$vectors %*% covewi %*% t(cove$vectors)
  }
  options(show.error.messages = TRUE)
  out <- list(inv = covinv, coll = coll)
  return(out)
}


#' Computes abundance estimates at specified parameter values using
#' Horvitz-Thompson-like estimator
#'
#' Computes abundance at specified values of parameters for numerical
#' computation of first derivative with respect to parameters in detection
#' function.  An internal function called by DeltaMethod which is invoked by
#' dht.se
#'
#' @param par detection function parameter values
#' @param model ddf model object
#' @param obs observations table
#' @param samples samples table
#' @param options list of options as specified in \code{\link{dht}}
#' @return vector of abundance estimates at values of parameters specified in
#'   par
#' @note Internal function; not intended to be called by user
#' @author Jeff Laake
#' @seealso \code{\link{dht}}, \code{\link{dht.se}}, \code{\link{DeltaMethod}}
#' @keywords utility
dht.deriv <- function(par,model,obs,samples,options=list()){
  # Functions Used:  predict, covered.region.dht, survey.region.dht
  
  #  Depending on model method store new parameter values
  #  uses coefficients instead of coef; needed for trial.fi and io.fi 
  model$par=par
  if(model$method!="ds"){
    if(model$method%in%c("io","trial","rem")){
      model$mr$mr$coefficients <- model$par[1:length(model$mr$mr$coefficients)]
      model$ds$par <-model$par[(length(model$mr$mr$coefficients)+1):length(par)]
    }else{
      model$mr$coefficients <- model$par
    }
  }
  
  # Get new predicted detection probabilities and put with observations
  # 1/1/05 jll; change made to fi methods. need integrate=FALSE
  # 1/24/06 jll; must have been on glue when I made the 1/1/05 change; clearly
  # not what I wanted.  You would only not integrate if you didn't want to
  # assume 1/w; otherwise, always want to integrate to get average detection
  # probability
  pdot <- predict(model,integrate=TRUE,compute=TRUE)
  
  if(!is.null(pdot$fitted)){
    pdot <- pdot$fitted
  }
  
  # probably should stop here if pdot$fitted is null? anyway...
  # code below would not work if there is a pdot column there already, so remove
  obs$pdot <- NULL
  
  obs <- merge(obs,data.frame(object=as.numeric(names(model$fitted)),pdot=pdot))
  
  # Compute covered region abundances by sample depending on value of group
  Nhat.by.sample <- covered.region.dht(obs,samples,options$group)
  
  # Scale up abundances to survey region
  Nhat.by.sample <- survey.region.dht(Nhat.by.sample, samples,
                                      model$meta.data$width*options$convert.units,
                                      model$meta.data$point)
  Nhat.by.region <- by(Nhat.by.sample$Nhat,Nhat.by.sample$Region.Label,sum)
  
  # Return vector of predicted abundances
  
  if(length(Nhat.by.region)>1){
    return(c(as.vector(Nhat.by.region),sum(as.vector(Nhat.by.region))))
  }else{
    return(as.vector(Nhat.by.region))
  }
}