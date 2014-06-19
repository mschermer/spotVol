#' Spot volatility estimation
#' 
#' The \code{spotvolatility} package offers several methods to estimate spot volatility
#' and its intraday seasonality, using high-frequency data.
#' 
#' @details
#' The following spot volatility estimation methods have been implemented:
#' 
#' @section Deterministic periodicity:
#' The spot volatility is decomposed into a a deterministic periodic factor f_{i}
#'  (identical for every day in the sample) and a daily factor s_{t} 
#'  (identical for all observations within a day). Both components are then
#'  estimated separately. For more details, see Taylor and Xu (1997) and 
#'  Andersen and Bollerslev (1997). The jump robust versions by Boudt et al. (2011)
#'  have also been implemented.
#' 
#' @section Stochastic periodicity:
#' This method by Beltratti and Morana (2001) assumes the periodicity factor to be stochastic.
#' The spot volatility estimation is split into four components: a random walk, an autoregressive
#' process, a stochastic cyclical process and a deterministic cyclical process. The model is 
#' estimated using a quasi-maximum likelihood method based on the Kalman Filter. The package
#' \code{FKF} is used to apply the Kalman filter.
#' 
#' @references
#' Andersen, T. G. and T. Bollerslev (1997). Intraday periodicity and volatility
#' persistence in financial markets. Journal of Empirical Finance 4, 115-158.
#' 
#' Beltratti, A. and C. Morana (2001). Deterministic and stochastic methods for estimation
#' of intraday seasonal components with high frequency data. Economic Notes 30, 205-234.
#' 
#' Boudt K., Croux C. and Laurent S. (2011). Robust estimation of intraweek periodicity
#' in volatility and jump detection. Journal of Empirical Finance 18, 353-367.
#' 
#' Kristensen, Dennis (2010). Nonparametric filtering of the realized spot volatility: 
#' A kernel-based approach. Econometric Theory 26, 60-93.
#' 
#' Taylor, S. J. and X. Xu (1997). The incremental volatility information in one million
#' foreign exchange quotations. Journal of Empirical Finance 4, 317-340.
#' 
#' @docType package
#' @name spotvolatility
#' @import highfrequency FKF chron timeDate
NULL

#' spotvol
#' 
#' @description
#'
#' The \code{spotvol} function offers several methods to estimate spot volatility
#' and its intraday seasonality, using high-frequency data. It returns an object of 
#' class \code{spotvol}, which can contain various outputs, depending on the method
#' used. See 'Details' for a description of each method. In any case, the output will
#' contain the spot volatility estimates. 
#' 
#' The input can consist of price data or return data, either tick by tick or sampled at
#' set intervals. The data will be converted to equispaced high-frequency returns \eqn{r_{t,i}}
#' (read: the \eqn{i}th return on day \eqn{t}). 
#' 
#' @details The following estimation methods can be specified in \code{method}:
#' 
#' \strong{Deterministic periodicity method (\code{"detper"})}
#' 
#' Parameters:
#' \tabular{ll}{
#' \code{dailyvol} \tab A string specifying the estimation method for the daily component 
#' \eqn{s_t}. Possible values are \code{"bipower", "rv", "medrv"}. Default = \code{"bipower"}. \cr
#' \code{periodicvol} \tab A string specifying the estimation method for the component of intraday
#' volatility, that depends in a deterministic way on the intraday time at which the return is observed.
#' Possible values are \code{"TML", "SD", "WSD", "OLS"}. See Boudt et al. (2011) for details. Default = \code{"TML"}.\cr
#' \code{P1} \tab A positive integer corresponding to the number of cosinus terms
#'  used in the flexible Fourier specification of the periodicity function, 
#'  see Andersen et al. (1997) for details. Default = 5. \cr
#' \code{P2} \tab Same as \code{P1}, but for the sinus terms. Default = 5.\cr
#' \code{dummies} \tab Boolean: in case it is \code{TRUE}, the parametric estimator of periodic standard
#' deviation specifies the periodicity function as the sum of dummy variables corresponding to each intraday period.
#' If it false, the parametric estimator uses the flexible Fourier specification. Default = \code{FALSE}.
#' }
#' Outputs (see 'Value' for a full description of each component):
#' \itemize{
#' \item{\code{spot}}
#' \item{\code{daily}}
#' \item{\code{periodic}}
#' }
#' The spot volatility is decomposed into a deterministic periodic factor \eqn{f_{i}} 
#' (identical for every day in the sample) and a daily factor \eqn{s_{t}} (identical for all
#' observations within a day). Both components are then estimated separately. For more
#' details, see Taylor and Xu (1997) and Andersen and Bollerslev (1997). The jump robust
#' versions by Boudt et al. (2011) have also been implemented.
#' 
#' \strong{Stochastic periodicity method (\code{"stochper"})}
#' 
#' Parameters:
#' \tabular{ll}{
#' \code{P1} \tab A positive integer corresponding to the number of cosinus terms
#'  used in the flexible Fourier specification of the periodicity function. 
#'  Default = 5. \cr
#' \code{P2} \tab Same as \code{P1}, but for the sinus terms. Default = 5.\cr
#' \code{init} \tab A named list of initial values to be used in the optimization routine (\code{"BFGS"}
#' in \code{optim}). Default = \code{list(sigma = 0.03, sigma_mu = 0.005, sigma_h = 0.005, sigma_k = 0.05, 
#' phi = 0.2, rho = 0.98, mu = c(2, -0.5), delta_c = rep(0, max(1,P1)), delta_s = rep(0, max(1,P2)))}.
#' See Beltratti & Morana (2001) for a definition of each parameter. \code{init} can contain any number of
#' these parameters. For parameters not specified in \code{init}, the default initial value will be used.\cr
#' \code{control} \tab A list of options to be passed down to \code{optim}.
#' }
#' Outputs (see 'Value' for a full description of each component):
#' \itemize{
#' \item{\code{spot}}
#' \item{\code{par}}
#' }
#' This method by Beltratti and Morana (2001) assumes the periodicity factor to be stochastic.
#' The spot volatility estimation is split into four components: a random walk, an autoregressive
#' process, a stochastic cyclical process and a deterministic cyclical process. The model is 
#' estimated using a quasi-maximum likelihood method based on the Kalman Filter. The package
#' \code{FKF} is used to apply the Kalman filter. In addition to the spot volatility estimates,
#' all parameter estimates are returned.
#' 
#' \strong{Nonparametric filtering (\code{"kernel"})}
#' 
#' Parameters:
#' \tabular{ll}{
#' \code{type} \tab String specifying the type of kernel to be used. Options include 
#' \code{"gaussian", "epanechnikov", "beta"}. Default = \code{"gaussian"}.\cr
#' \code{h} \tab Scalar or vector specifying bandwidth(s) to be used in kernel. If \code{h}
#' is a scalar, it will be assumed equal throughout the sample. If it is a vector, it should contain
#' bandwidths for each day. If left empty, it will be estimated. Default = \code{NULL}. \cr
#' \code{cv} \tab Boolean indicating whether cross-validation should be used to estimate
#' the bandwidth \code{h}. \code{cv} is overridden if \code{h} has already been specified 
#' by the user. Default = \code{TRUE}.\cr
#' \code{control} \tab A list of options to be passed down to \code{optim}.
#' }
#' Outputs (see 'Value' for a full description of each component):
#' \itemize{
#' \item{\code{spot}}
#' \item{\code{par}}
#' }
#' This method by Kristensen (2010) filters the spot volatility in a nonparametric way by applying
#' kernel weights to the standard realized volatility estimator. Different kernels and bandwidths can 
#' be used to focus on specific characteristics of the volatility process. When using this method, in
#' addition to the spot volatility estimates, all used values of the bandwidth \eqn{h} are returned.
#' 
#' @param data \code{xts} object, containing a price or return series. If the data consists 
#' of returns, set \code{makeReturns} to \code{FALSE}.
#' @param makeReturns boolean, if \code{TRUE} the function will calculate returns from the price
#' series \code{data}. If \code{data} already consists of returns, set to \code{FALSE}.
#' @param method specifies which method will be used to estimate the spot volatility. Options
#' include \code{"detper"} and \code{"stochper"}. See Details.
#' @param on string indicating the time scale in which \code{k} is expressed.
#'  Possible values are: \code{"secs", "seconds", "mins", "minutes", "hours"}.
#' @param k positive integer, indicating the number of periods to aggregate over.
#'  E.g. to aggregate an \code{xts} object to the 5 minute frequency, set \code{k = 5} and
#'  \code{on = "minutes"}.
#' @param marketopen the market opening time. By default, \code{marketopen = "09:30:00"}.
#' @param marketclose the market closing time. By default, \code{marketclose = "16:00:00"}.
#' @param ... method-specific parameters (see 'Details').
#' 
#' @return
#' A \code{spotvol} object, which is a list containing one or more of the following outputs,
#' depending on the method used:
#' 
#' \code{spot}
#' 
#' An \code{xts} object containing spot volatility estimates \eqn{\sigma_{t,i}}, reported for each interval
#' \eqn{i} between \code{marketopen} and \code{marketclose} for every day \eqn{t} in \code{data}. The length
#' of the intervals is specifiedby \code{k} and \code{on}. Methods that provide this output: All.
#' 
#' \code{daily}
#' 
#' An \code{xts} object containing estimates of the daily volatility levels for each day \eqn{t} in \code{data}, 
#' if the used method decomposed spot volatility into a daily and an intraday component. Methods that provide this
#' output: \code{"detper"}.
#' 
#' \code{periodic}
#' 
#' Estimates of the intraday periodicity factor for each day interval \eqn{i} between \code{marketopen}
#' and \code{marketclose}, if the spot volatility was decomposed into a daily and an intraday component.
#' Methods that provide this output: \code{"detper"}.
#' 
#' \code{par}
#' 
#' A named list containing parameter estimates, for methods that estimate one or more parameters. 
#' Methods that provide this output: \code{"stochper", "kernel"}.
#' 
#' 
#' @export
#' @examples
#' data(sample_prices_5min)
#' vol1 <- spotvol(sample_prices_5min)
#' init = list(sigma = 0.03, sigma_mu = 0.005, sigma_h = 0.007,
#'    sigma_k = 0.06, phi = 0.194, rho = 0.986, mu = c(1.87,-0.42),
#'    delta_c = c(0.25, -0.05, -0.2, 0.13, 0.02), delta_s = c(-1.2, 0.11, 0.26, -0.03, 0.08))
#' # next method will take around 110 iterations
#' vol2 <- spotvol(sample_prices_5min, method = "stochper", init = init)
#' vol3 <- spotvol(sample_prices_5min, method = "kernel", h = 0.05)
#' plot(vol1$spot, type="l")
#' lines(vol2$spot, col="red")
#' lines(vol3$spot, col="blue")
#' legend("topright", c("detper", "stochper", "kernel"), col = c("black", "red", "blue"), lty=1, inset=0.02, bg="white")
#' 
#' @references
#' Andersen, T. G. and T. Bollerslev (1997). Intraday periodicity and volatility
#' persistence in financial markets. Journal of Empirical Finance 4, 115-158.
#' 
#' Beltratti, A. and C. Morana (2001). Deterministic and stochastic methods for estimation
#' of intraday seasonal components with high frequency data. Economic Notes 30, 205-234.
#' 
#' Boudt K., Croux C. and Laurent S. (2011). Robust estimation of intraweek periodicity
#' in volatility and jump detection. Journal of Empirical Finance 18, 353-367.
#' 
#' Kristensen, Dennis (2010). Nonparametric filtering of the realized spot volatility: 
#' A kernel-based approach. Econometric Theory 26, 60-93.
#' 
#' Taylor, S. J. and X. Xu (1997). The incremental volatility information in one million
#' foreign exchange quotations. Journal of Empirical Finance 4, 317-340.
spotvol =  function(data, makeReturns = TRUE, method = "detper", on = "minutes", k = 5,
                    marketopen = "09:30:00", marketclose = "16:00:00", ...)  
{
  Sys.setenv(TZ="GMT") # only for convenience during testing
  dates = unique(format(time(data), "%Y-%m-%d"))
  cDays = length(dates)
  rdata = mR = c()
  if(on=="minutes"){
    intraday = seq(from=times(marketopen), to=times(marketclose), by=times(paste("00:0",k,":00",sep=""))) 
  }
  if(as.character(tail(intraday,1))!=marketclose) intraday=c(intraday,marketclose)
  if(makeReturns) intraday = intraday[2:length(intraday)]
  for (d in 1:cDays) {
    datad = data[as.character(dates[d])]
    datad = aggregatePrice(datad, on = on, k = k , marketopen = marketopen, marketclose = marketclose)
    z = xts( rep(1,length(intraday)) , order.by = timeDate( paste(dates[d],as.character(intraday),sep="") , format = "%Y-%m-%d %H:%M:%S"))
    datad = merge.xts( z , datad )$datad
    datad = na.locf(datad)
    if (makeReturns)
    {
      rdatad = makeReturns(datad)
      rdatad = rdatad[time(rdatad) > min(time(rdatad))]
      rdata = rbind(rdata, rdatad)
      mR = rbind(mR, as.numeric(rdatad))
    }
    else 
    {
      rdata = rbind(rdata, datad)
      mR = rbind(mR, as.numeric(datad))
    } 
  }
  mR[is.na(mR)]=0
  options <- list(...)
  out = switch(method, 
           detper = detper(mR, options = options), 
           stochper = stochper(mR, options = options),
           kernel = kernelestim(mR, options = options))  
  if("periodic" %in% names(out)) names(out$periodic) <- as.character(intraday)
  
  return(out)
}

# Deterministic periodicity model
# 
# Modified spotVol function from highfrequency package
detper = function(mR, options = list()) 
{
  # default options, replace if user-specified
  op <- list(dailyvol = "bipower", periodicvol = "TML", dummies = FALSE, P1 = 5, P2 = 5)
  op[names(options)] <- options 
  
  cDays = nrow(mR)
  M = ncol(mR)
  if (cDays == 1) { 
    mR = as.numeric(mR)
    estimdailyvol = switch(op$dailyvol, 
                           bipower = rBPCov(mR), 
                           medrv = medRV(mR), rv = RV(mR))
  }else {
    estimdailyvol = switch(op$dailyvol, 
                           bipower = apply(mR, 1, "rBPCov"),
                           medrv = apply(mR, 1, "medRV"), rv = apply(mR, 1, "RV"))
  }  
  if (cDays <= 50) {
    print("Periodicity estimation requires at least 50 observations. Periodic component set to unity")
    estimperiodicvol = rep(1, M)
  }
  else {
    mstdR = mR/sqrt(estimdailyvol * (1/M))
    selection = c(1:M)[ (nrow(mR)-apply(mR,2,'countzeroes')) >=20] 
    # preferably no na is between
    selection = c( min(selection) : max(selection) )
    mstdR = mstdR[,selection]
    estimperiodicvol_temp = diurnal(stddata = mstdR, method = op$periodicvol, 
                                    dummies = op$dummies, P1 = op$P1, P2 = op$P2)[[1]]
    estimperiodicvol = rep(1,M)
    estimperiodicvol[selection] = estimperiodicvol_temp
    mfilteredR = mR/matrix(rep(estimperiodicvol, cDays), 
                         byrow = T, nrow = cDays)
    estimdailyvol = switch(op$dailyvol, bipower = apply(mfilteredR, 1, "rBPCov"),
                            medrv = apply(mfilteredR, 1, "medRV"), 
                            rv = apply(mfilteredR, 1, "RV"))
    out <- list(spot = rep(sqrt(estimdailyvol * (1/M)), each = M) * rep(estimperiodicvol, cDays),
                daily = estimdailyvol,
                periodic = estimperiodicvol)
    class(out) <- "spotvol"
    return(out)
  }
}

# Stochastic periodicity model
# 
# This function estimates the spot volatility by using the stochastic periodcity model of Beltratti & Morana (2001)
stochper =  function(mR, options = list()) 
{
  # default options, replace if user-specified
  op <- list(init = list(), P1 = 5, P2 = 5, control = list(trace=1, maxit=500))
  op[names(options)] <- options 
  
  N = ncol(mR)
  days = nrow(mR)
  
  logr2 = log(mR^2)
  rvector = as.vector(t(logr2)) 
  lambda = (2*pi)/N;
  
  # default starting values of parameters
  sp <- list(sigma = 0.03,
    sigma_mu = 0.005,
    sigma_h = 0.005,
    sigma_k = 0.05,
    phi = 0.2,
    rho = 0.98,
    mu = c(2, -0.5),
    delta_c = rep(0, max(1,op$P1)),
    delta_s = rep(0, max(1,op$P2)))
  
  # replace if user has specified different values
  sp[names(op$init)] <- op$init
  
  # check input
  for (i in c("sigma", "sigma_mu", "sigma_h", "sigma_k", "phi", "rho"))
  {
      if (sapply(sp, length)[i] != 1) stop(paste(i, " must be a scalar"))  
  }
  if (length(sp$mu) != 2) stop("mu must have length 2")
  if (length(sp$delta_c) != op$P1 & op$P1 > 0) stop("delta_c must have length equal to P1")
  if (length(sp$delta_s) != op$P2 & op$P2 > 0) stop("delta_s must have length equal to P2")
  if (length(sp$delta_c) < 1) stop("delta_c must at least have length 1")
  if (length(sp$delta_s) < 1) stop("delta_s must at least have length 1")
  
  # transform parameters to allow for unrestricted optimization (domain -Inf to Inf)
  par_t <- c(sigma = log(sp$sigma), sigma_mu = log(sp$sigma_mu), sigma_h = log(sp$sigma_h),
               sigma_k = log(sp$sigma_k), phi = log(sp$phi/(1-sp$phi)), rho = log(sp$rho/(1-sp$rho)),
               mu = sp$mu, delta_c = sp$delta_c, delta_s = sp$delta_s)   

  opt <- optim(par_t, loglikBM, yt = rvector, N = N, days = days, P1 = op$P1, P2 = op$P2, method="BFGS", control = op$control)
  
  # recreate model to obtain volatility estimates
  ss <- ssmodel(opt$par, days, N, P1 = op$P1, P2 = op$P2)
  kf <- fkf(a0 = ss$a0, P0 = ss$P0, dt = ss$dt, ct = ss$ct, Tt = ss$Tt, Zt = ss$Zt, 
            HHt = ss$HHt, GGt = ss$GGt, yt = matrix(rvector, ncol = length(rvector)))
  sigmahat <- as.vector(exp((ss$Zt%*%kf$at[,1:(N*days)] + ss$ct + 1.27)/2))
  
  # transform parameter estimates back
  estimates <- c(exp(opt$par["sigma"]), exp(opt$par["sigma_mu"]), exp(opt$par["sigma_h"]), exp(opt$par["sigma_k"]),
                 exp(opt$par["phi"])/(1+exp(opt$par["phi"])), exp(opt$par["rho"])/(1+exp(opt$par["rho"])), opt$par[-(1:6)])

  out <- list(spot = sigmahat,
              par = estimates)
  class(out) <- "spotvol"
  return(out)
}

# Calculate log likelihood using Kalman Filter
# 
# This function returns the average log likehood value of the stochastic periodicity model, given the input parameters.
loglikBM <- function(par_t, yt, days, N = 288, P1 = 5, P2 = 5)
{
  ss <- ssmodel(par_t, days, N, P1 = P1, P2 = P2)
  yt = matrix(yt, ncol = length(yt))
  kf <- fkf(a0 = ss$a0, P0 = ss$P0, dt = ss$dt, ct = ss$ct, Tt = ss$Tt, Zt = ss$Zt, HHt = ss$HHt, GGt = ss$GGt, yt = yt)
  return(-kf$logLik/length(yt))
}

# Generate state space model
# 
# This function creates the state space matrices from the input parameters.
# The output is in the format used by the FKF package.
ssmodel <- function(par_t, days, N = 288, P1 = 5, P2 = 5)
{
  par <- c(exp(par_t["sigma"]), exp(par_t["sigma_mu"]), exp(par_t["sigma_h"]), exp(par_t["sigma_k"]), 
           exp(par_t["phi"])/(1+exp(par_t["phi"])), exp(par_t["rho"])/(1+exp(par_t["rho"])), par_t[-(1:6)])
  lambda=(2*pi)/288
  a0 <- c(0, 0, par["delta_c1"], par["delta_s1"])
  if (P1 == 0) a0[3] <- par["delta_c"]
  if (P2 == 0) a0[4] <- par["delta_s"]   
  m <- length(a0)
  P0 <- Tt <- Ht <- matrix(0, m, m)
  diag(Tt) <- c(1, par["phi"], rep(par["rho"]*cos(lambda), 2))
  Tt[3,4] <- par["rho"]*sin(lambda)
  Tt[4,3] <- par["rho"]*-sin(lambda)
  Zt <- matrix(c(1, 1, 1, 0), ncol=m)
  Gt <- sqrt(0.5*pi^2)
  GGt <- Gt %*% t(Gt)
  diag(Ht) <- c(par["sigma_mu"], par["sigma_h"], rep(par["sigma_k"], 2))
  HHt <- Ht %*% t(Ht)
  dt <- matrix(0, nrow=m)
  ct <- log(par["sigma"]^2) - 1.270363
  
  # calculate deterministic part c2, add to ct
  n <- 1:N
  M1 <- (2*n)/(N+1)
  M2 <- (6*n^2)/((N+1)*(N+2))
  c2 <- par["mu1"]*M1 + par["mu2"]*M2
  if (P1 > 1)
  {
    for (k in 2:P1)
    {
      c2 <- c2 + par[paste("delta_c", k, sep="")]*cos(k*lambda*n) 
    }  
  }
  if (P2 > 1)
  {
    for (p in 2:P2)
    {
      c2 <- c2 + par[paste("delta_s", p, sep="")]*sin(p*lambda*n)
    }  
  }
  ct <- matrix(ct + c2, ncol = N*days)
  
  return(list(a0 = a0, P0 = P0, Tt = Tt, Zt = Zt, GGt = GGt, HHt = HHt, dt = dt, ct = ct))
}

# Kernel estimation method
# 
# See Kristensen (2010)
kernelestim <- function(mR, options = list())
{
  # default options, replace if user-specified
  op <- list(type = "gaussian", h = NULL, cv = TRUE, control = list())
  op[names(options)] <- options 
  
  D = nrow(mR)
  N = ncol(mR)
  t <- (1:N)/N
  if (is.null(op$h)) h <- numeric(D)
  else h <- rep(op$h, length.out = D)
  
  sigma2hat <- matrix(NA, nrow = D, ncol = N)
  for(d in 1:D)
  {
    if (is.null(op$h))
    {
      cat(paste("Estimating optimal bandwidth for day", d, "of", D, "...\n"))
      h[d] <- estbandwidth(mR[d, ], type = op$type, cv = op$cv, control = op$control)
    }
    for(n in 1:N)
    {
      if (op$type == "beta") 
      {
        K <- sapply(t, 'kernelk', type = "beta", b = h[d], y = t[n])
      }
      else
      {
        K <- (sapply((t - t[n])/h[d], 'kernelk', type = op$type))/h[d] 
      }
      K <- K/sum(K)
      sigma2hat[d, n] <- K %*% (mR[d, ]^2)
    }
  }
  out = list(spot = as.vector(t(sqrt(sigma2hat))), par = list(h = h))
  class(out) <- "spotvol"
  return(out)
}

# calculate values of certain kernels
# arguments b and n only needed for type == "beta"
kernelk <- function(x, type = "gaussian", b = 1, y = 1)
{
  if (type == "gaussian") return(dnorm(x))  
  if (type == "epanechnikov")
  {
    if (abs(x) > 1) return(0)
    else return((3/4)*(1-x)^2)
  }
  if (type == "beta") return(dbeta(x, y/b + 1, (1-y)/b + 1))
}

# estimate optimal bandwidth paramater h
# by default, this is done through crossvalidation (cv)
# else the formula for h_opt in Kristensen(2010) is approximated
estbandwidth <- function(x, type = "gaussian", cv = TRUE, control = list())
{
  if (cv)
  {
    opt <- optim(log(0.05), ISE, method = "BFGS", control = control, x = x, type = type)
    h = exp(opt$par)
  }
  else
  {
    N = length(x)
    Delta = 1/N
    quarticity = (Delta/3)*sum(x^4)
    if (type == "gaussian")
    {
      q = 2
      RK = 1/(2*pi)
      kq2 = 1
    }
    if (type == "epanechnikov")
    {
      q = 2
      RK = 3/5
      kq2 = 1/5
    }
    h = (((quarticity*RK)/(q*kq2))^(1/(2*q+1)))*(N^(-1/(2*q+1)))
  }
  return(h)
}

ISE <- function(h, x, type = "gaussian")
{
  h = exp(h)
  N = length(x)
  t = (1:N)/N
  sigma2hat <- rep(NA, N)
  for(n in 1:N)
  {
    if (type == "beta") 
    {
      K <- sapply(t, 'kernelk', type = "beta", b = h, y = t[n])
    }
    else
    {
      K <- (sapply((t - t[n])/h, 'kernelk', type = type))/h
    }
    K[n] = 0
    K <- K/sum(K)
    sigma2hat[n] <- K %*% (x^2) 
  }    
  tl = 5
  tu = N-5
  ISE = sum(abs((x[tl:tu]^2) - sigma2hat[tl:tu]))
  return(ISE)
}


### auxiliary internal functions copied from highfrequency package

countzeroes = function( series )
{
  return( sum( 1*(series==0) ) )
}

HRweight = function( d,k){
  # Hard rejection weight function
  w = 1*(d<=k); return(w)
}

shorthscale = function( data )
{
  sorteddata = sort(data);
  n = length(data);
  h = floor(n/2)+1;
  M = matrix( rep(0,2*(n-h+1) ) , nrow= 2 );
  for( i in 1:(n-h+1) ){
    M[,i] = c( sorteddata[ i ], sorteddata[ i+h-1 ] )
  }
  return( 0.7413*min( M[2,]-M[1,] ) );
}

diurnal = 
  function (stddata, method = "TML", dummies = F, P1 = 6, P2 = 4) 
  {
    cDays = dim(stddata)[1]
    intraT = dim(stddata)[2]
    meannozero = function(series) {
      return(mean(series[series != 0]))
    }
    shorthscalenozero = function(series) {
      return(shorthscale(series[series != 0]))
    }
    WSDnozero = function(weights, series) {
      out = sum((weights * series^2)[series != 0])/sum(weights[series != 
                                                                 0])
      return(sqrt(1.081 * out))
    }
    if (method == "SD" | method == "OLS") {
      seas = sqrt(apply(stddata^2, 2, "meannozero"))
    }
    if (method == "WSD" | method == "TML") {
      seas = apply(stddata, 2, "shorthscalenozero")
      shorthseas = seas/sqrt(mean(seas^2))
      shorthseas[shorthseas == 0] = 1
      weights = matrix(HRweight(as.vector(t(stddata^2)/rep(shorthseas, 
                                                           cDays)^2), qchisq(0.99, df = 1)), ncol = dim(stddata)[2], 
                       byrow = T)
      for (c in 1:intraT) {
        seas[c] = WSDnozero(weights[, c], stddata[, c])
      }
    }
    seas = na.locf(seas,na.rm=F) #do not remove leading NA
    seas = na.locf(seas,fromLast=T)
    seas = seas/sqrt(mean(seas^2))
    if (method == "OLS" | method == "TML") {
      c = center()
      vstddata = as.vector(stddata)
      nobs = length(vstddata)
      vi = rep(c(1:intraT), each = cDays)
      if (method == "TML") {
        if( length(vstddata)!= length(seas)*cDays ){ print(length(vstddata)); print(length(seas)); print(cDays)}
        firststepresids = log(abs(vstddata)) - c - log(rep(seas, 
                                                           each = cDays))
      }
      X = c()
      if (!dummies) {
        if (P1 > 0) {
          for (j in 1:P1) {
            X = cbind(X, cos(2 * pi * j * vi/intraT))
          }
        }
        M1 = (intraT + 1)/2
        M2 = (2 * intraT^2 + 3 * intraT + 1)/6
        ADD = (vi/M1)
        X = cbind(X, ADD)
        ADD = (vi^2/M2)
        X = cbind(X, ADD)
        if (P2 > 0) {
          ADD = c()
          for (j in 1:P2) {
            ADD = cbind(ADD, sin(2 * pi * j * vi/intraT))
          }
        }
        X = cbind(X, ADD)
        opening = vi - 0
        stdopening = (vi - 0)/80
        almond1_opening = (1 - (stdopening)^3)
        almond2_opening = (1 - (stdopening)^2) * (opening)
        almond3_opening = (1 - (stdopening)) * (opening^2)
        X = cbind(X, almond1_opening, almond2_opening, almond3_opening)
        closing = max(vi) - vi
        stdclosing = (max(vi) - vi)/max(vi)
        almond1_closing = (1 - (stdclosing)^3)
        almond2_closing = (1 - (stdclosing)^2) * (closing)
        almond3_closing = (1 - (stdclosing)) * (closing^2)
        X = cbind(X, almond1_closing, almond2_closing, almond3_closing)
      }
      else {
        for (d in 1:intraT) {
          dummy = rep(0, intraT)
          dummy[d] = 1
          dummy = rep(dummy, each = cDays)
          X = cbind(X, dummy)
        }
      }
      selection = c(1:nobs)[vstddata != 0]
      vstddata = vstddata[selection]
      X = X[selection, ]
      if (method == "TML") {
        firststepresids = firststepresids[selection]
      }
      vy = matrix(log(abs(vstddata)), ncol = 1) - c
      if (method == "OLS") {
        Z = try(solve(t(X) %*% X), silent = T)
        if (inherits(Z, "try-error")) {
          print("X'X is not invertible. Switch to TML")
        }
        else {
          theta = solve(t(X) %*% X) %*% t(X) %*% vy
          rm(X)
          rm(vy)
        }
      }
      if (method == "TML") {
        inittheta = rep(0, dim(X)[2])
        l = -2.272
        u = 1.6675
        nonoutliers = c(1:length(vy))[(firststepresids > 
                                         l) & (firststepresids < u)]
        truncvy = vy[nonoutliers]
        rm(vy)
        truncX = X[nonoutliers, ]
        rm(X)
        negtruncLLH = function(theta) {
          res = truncvy - truncX %*% matrix(theta, ncol = 1)
          return(mean(-res - c + exp(2 * (res + c))/2))
        }
        grnegtruncLLH = function(theta) {
          res = truncvy - truncX %*% matrix(theta, ncol = 1)
          dres = -truncX
          return(apply(-dres + as.vector(exp(2 * (res + 
                                                    c))) * dres, 2, "mean"))
        }
        est = optim(par = inittheta, fn = negtruncLLH, gr = grnegtruncLLH, 
                    method = "BFGS")
        theta = est$par
        rm(truncX)
        rm(truncvy)
      }
# disable plot for now      
#       plot(seas, main = "Non-parametric and parametric periodicity estimates", 
#            xlab = "intraday period", type = "l", lty = 3)
#       legend("topright", c("Parametric", "Non-parametric"), cex = 1.1,
#              lty = c(1,3), lwd = 1, bty = "n")
        seas = highfrequency:::diurnalfit(theta = theta, P1 = P1, P2 = P2, intraT = intraT, 
                                         dummies = dummies)
#       lines(seas, lty = 1)
      return(list(seas, theta))
    }
    else {
      return(list(seas))
    }
  }

diurnalfit = function( theta , P1 , P2 , intraT , dummies=F )
{
  vi = c(1:intraT) ;  
  M1 = (intraT+1)/2 ; M2 = (2*intraT^2 + 3*intraT + 1)/6;
  
  # Regressors that do not depend on Day of Week:
  X = c()
  if(!dummies){
    if ( P1 > 0 ){ for( j in 1:P1 ){ X = cbind( X , cos(2*pi*j*vi/intraT) )   }  } 
    
    ADD = (vi/M1 ) ; X = cbind(X,ADD);
    ADD = (vi^2/M2); X = cbind(X,ADD);
    if ( P2 > 0 ){ ADD= c(); for( j in 1:P2 ){  ADD = cbind( ADD , sin(2*pi*j*vi/intraT)  ) }}; X = cbind( X , ADD ) ; 
    
    #openingeffect
    opening = vi-0 ; stdopening = (vi-0)/80 ;
    almond1_opening   = ( 1 - (stdopening)^3 );
    almond2_opening   = ( 1 - (stdopening)^2 )*( opening);
    almond3_opening   = ( 1 - (stdopening)   )*( opening^2);   
    X = cbind(  X, almond1_opening , almond2_opening , almond3_opening   )  ;
    
    #closing effect
    closing = max(vi)-vi ; stdclosing = (max(vi)-vi)/max(vi) ;
    almond1_closing   = ( 1 - (stdclosing)^3 );
    almond2_closing   = ( 1 - (stdclosing)^2 )*( closing);
    almond3_closing   = ( 1 - (stdclosing)   )*( closing^2);   
    X = cbind(  X, almond1_closing , almond2_closing , almond3_closing   )  ;
    
  }else{
    for( d in 1:intraT){
      dummy = rep(0,intraT); dummy[d]=1; 
      X = cbind(X,dummy); 
    }
  }
  # Compute fit
  seas = exp( X%*%matrix(theta,ncol=1) );
  seas = seas/sqrt(mean( seas^2) )    
  return( seas )          
}

center = function()
{
  g=function(y){ return( sqrt(2/pi)*exp(y-exp(2*y)/2)  )}
  f=function(y){ return( y*g(y)    )  }
  return( integrate(f,-Inf,Inf)$value )
}

#' Sample prices data
#' 
#' 'sample_real5minprices' from highfrequency package
#' 
#' @docType data
#' @name sample_prices_5min
#' @usage data(sample_prices_5min)
#' @format xts
NULL

#' Sample returns data
#' 
#' EUR/USD returns from January to September 2004
#' 
#' @docType data
#' @name sample_returns_5min
#' @usage data(sample_returns_5min)
#' @format xts
NULL