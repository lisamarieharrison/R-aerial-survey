#functions from mrds package for fiddling with gamma detection function apex

myDeltaMethod <- function (par, fct, vcov, delta, ...) 
{
  theta <- function(par) fct(par, ...)
  savepar <- par
  value1 <- theta(par)
  partial <- matrix(0, nrow = length(par), ncol = length(value1))
  for (i in 1:length(par)) {
    if (savepar[i] != 0) {
      deltap <- delta * savepar[i]
    }
    else {
      deltap <- delta
    }
    par <- savepar
    par[i] <- savepar[i] + deltap
    value1 <- theta(par)
    par <- savepar
    par[i] <- savepar[i] - deltap
    value2 <- theta(par)
    partial[i, ] <- (value1 - value2)/(2 * deltap)
  }
  variance <- t(partial) %*% vcov %*% partial
  print(fct)
  return(list(variance = variance, partial = partial))
}

fct <- function (par, model = NULL, group = TRUE) 
{
  if (is.null(model)) {
    model <- par
    par <- NULL
  }
  result <- switch(model$method, ds = NCovered.ds(par = par, 
                                                  model = model, group = group), io = NCovered.io(par = par, 
                                                                                                  model = model, group = group), io.fi = NCovered.io.fi(par = par, 
                                                                                                                                                        model = model, group = group), trial = NCovered.trial(par = par, 
                                                                                                                                                                                                              model = model, group = group), trial.fi = NCovered.trial.fi(par = par, 
                                                                                                                                                                                                                                                                          model = model, group = group), rem = NCovered.rem(par = par, 
                                                                                                                                                                                                                                                                                                                            model = model, group = group), rem.fi = NCovered.rem.fi(par = par, 
                                                                                                                                                                                                                                                                                                                                                                                    model = model, group = group))
  return(result)
}

NCovered.trial <- function (par, model, group = TRUE, ...) 
{
  if (!is.null(par)) {
    model$mr$mr$coefficients <- par[1:length(model$mr$mr$coefficients)]
    model$ds$par <- par[(length(model$mr$mr$coefficients) + 
                           1):length(par)]
    model$ds$ds$aux$ddfobj <- assign.par(model$ds$ds$aux$ddfobj, 
                                         model$ds$par)
    fitted <- predict(model, compute = TRUE, integrate = TRUE)$fitted
  }
  else {
    fitted <- model$fitted
  }
  if (!group) {
    size <- model$data$size[model$data$observer == 1 & model$data$object %in% 
                              as.numeric(names(model$fitted))]
    Nhat <- sum(compute.Nht(fitted, FALSE, size))
  }
  else {
    Nhat <- sum(compute.Nht(fitted, TRUE, size = NULL))
  }
  return(Nhat)
}

assign.par <- function (ddfobj, fpar) 
{
  index <- 1
  if (!is.null(ddfobj$shape)) {
    ddfobj$shape$parameters <- fpar[1:ncol(ddfobj$shape$dm)]
    index <- index + ncol(ddfobj$shape$dm)
  }
  if (!is.null(ddfobj$scale)) {
    ddfobj$scale$parameters <- fpar[index:(index + ncol(ddfobj$scale$dm) - 
                                             1)]
    index <- index + ncol(ddfobj$scale$dm)
  }
  if (!is.null(ddfobj$adjustment)) {
    ddfobj$adjustment$parameters <- fpar[index:length(fpar)]
  }
  return(ddfobj)
}

calc.se.Np <- function (model, avgp, n, average.p) 
{
  se.obj <- list()
  vcov <- solvecov(model$hessian)$inv
  Nhatvar.list <- DeltaMethod(model$par, NCovered, vcov, 0.001, 
                              model = model, group = TRUE)
  Nhatvar <- Nhatvar.list$variance + sum((1 - model$fitted)/model$fitted^2)
  cvN <- sqrt(Nhatvar)/model$Nhat
  var.pbar.list <- prob.se(model, avgp, vcov)
  covar <- t(Nhatvar.list$partial) %*% vcov %*% var.pbar.list$partial + 
    var.pbar.list$covar
  var.pbar <- average.p^2 * (cvN^2 + var.pbar.list$var/n^2 - 
                               2 * covar/(n * model$Nhat))
  se.obj$Nhat.se <- sqrt(Nhatvar)
  se.obj$average.p.se <- sqrt(var.pbar)
  se.obj$Nhatvar.list <- Nhatvar.list
  se.obj$vcov <- vcov
  return(se.obj)
}


NCovered.trial.fi <- function (par = NULL, model, group = TRUE, ...) 
{
  if (!is.null(par)) {
    model$mr$coefficients <- par
    fitted <- predict(model, compute = TRUE, integrate = TRUE)$fitted
  }
  else {
    fitted <- model$fitted
  }
  if (!group) {
    size <- model$data$size[model$data$observer == 1 & model$data$object %in% 
                              as.numeric(names(model$fitted))]
    Nhat <- sum(compute.Nht(fitted, FALSE, size))
  }
  else {
    Nhat <- sum(compute.Nht(fitted, TRUE, size = NULL))
  }
  return(Nhat)
}

prob.se <- function (model, fct, vcov, observer = NULL, fittedmodel = NULL) 
{
  if (is.null(fittedmodel)) {
    vc1.list <- DeltaMethod(model$par, prob.deriv, vcov, 
                            1e-04, model = model, parfct = fct, observer = observer, 
                            fittedmodel = NULL)
  }
  else {
    vc1.list <- DeltaMethod(fittedmodel$par, prob.deriv, 
                            vcov, 1e-04, model = model, parfct = fct, observer = observer, 
                            fittedmodel = fittedmodel)
  }
  vc1 <- vc1.list$variance
  if (!is.null(observer)) {
    newdat <- model$mr$data
    newdat$distance <- rep(0, length(newdat$distance))
    newdat$offsetvalue <- 0
    pred.at0 <- predict(model, newdat)
  }
  else if (!is.null(fittedmodel)) {
    if (class(fittedmodel)[1] != "rem") {
      newdat <- model$data[model$data$observer == 1 & model$data$detected == 
                             1, ]
    }
    else {
      newdat <- model$data
    }
    newdat <- newdat[newdat$distance <= model$meta.data$width & 
                       newdat$distance >= model$meta.data$left, ]
    newdat$distance <- rep(0, length(newdat$distance))
    newdat$offsetvalue <- 0
    pred.at0 <- predict(model, newdat)$fitted
  }
  if (is.null(fittedmodel)) {
    pdot <- model$fitted
    if (is.null(observer)) {
      vc2 <- sum(fct(model, pdot)^2 * (1 - pdot)/pdot^2)
      covar <- sum(fct(model, pdot) * (1 - pdot)/pdot^2)
    }
    else {
      vc2 <- sum(fct(model, pred.at0, observer)^2 * (1 - 
                                                       pdot)/pdot^2)
      covar <- sum(fct(model, pred.at0, observer) * (1 - 
                                                       pdot)/pdot^2)
    }
  }
  else {
    pdot <- fittedmodel$fitted
    vc2 <- sum(fct(model, pred.at0, observer)^2 * (1 - pdot)/pdot^2)
    covar <- sum(fct(model, pred.at0, observer) * (1 - pdot)/pdot^2)
  }
  return(list(var = vc1 + vc2, partial = vc1.list$partial, 
              covar = covar))
}


prob.deriv <- function (par, model, parfct, observer = NULL, fittedmodel = NULL) 
{
  set.par <- function(model, par) {
    model$par <- par
    if (model$method != "ds") {
      if (model$method == "io" | model$method == "trial" | 
          model$method == "rem") {
        model$mr$mr$coefficients <- model$par[1:length(model$mr$mr$coefficients)]
        model$ds$par <- model$par[(length(model$mr$mr$coefficients) + 
                                     1):length(par)]
        model$ds$ds$aux$ddfobj <- assign.par(model$ds$ds$aux$ddfobj, 
                                             model$ds$par)
      }
      else {
        model$mr$coefficients <- model$par
      }
    }
    return(model)
  }
  model <- set.par(model, par)
  if (!is.null(observer)) {
    newdat <- model$mr$data
    newdat$distance <- rep(0, length(newdat$distance))
    newdat$offsetvalue <- 0
    pred.at0 <- predict(model, newdat)
  }
  else {
    if (!is.null(fittedmodel)) {
      if (class(fittedmodel)[1] != "rem") {
        newdat <- model$data[model$data$observer == 1 & 
                               model$data$detected == 1, ]
      }
      else {
        newdat <- model$data
      }
      newdat <- newdat[newdat$distance <= model$meta.data$width & 
                         newdat$distance >= model$meta.data$left, ]
      newdat$distance <- rep(0, length(newdat$distance))
      newdat$offsetvalue <- 0
      pred.at0 <- predict(model, newdat)$fitted
    }
  }
  if (is.null(fittedmodel)) {
    pdot <- predict(model, compute = TRUE)$fitted
    if (is.null(observer)) {
      return(sum(parfct(model, pdot)/pdot))
    }
    else {
      return(sum(parfct(model, pred.at0, observer)/pdot))
    }
  }
  else {
    fittedmodel <- set.par(fittedmodel, par)
    pdot <- predict(fittedmodel, compute = TRUE)$fitted
    if (is.null(observer)) {
      return(sum(parfct(model, pred.at0)/pdot))
    }
    else {
      return(sum(parfct(model, pred.at0, observer)/pdot))
    }
  }
}