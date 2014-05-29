#' spotvol
#'
#' Estimating Spot Volatility
#' 
#' @export
spotvol =  function(pdata, dailyvol = "bipower", periodicvol = "TML", on = "minutes", 
                    k = 5, dummies = FALSE, P1 = 4, P2 = 2,  marketopen = "09:30:00", 
                    marketclose = "16:00:00") 
{
  require(chron);
  require(timeDate);
  require(highfrequency);
  dates = unique(format(time(pdata), "%Y-%m-%d"))
  cDays = length(dates)
  rdata = mR = c()
  if(on=="minutes"){
    intraday = seq(from=times(marketopen), to=times(marketclose), by=times(paste("00:0",k,":00",sep=""))) 
  }
  if(tail(intraday,1)!=marketclose){intraday=c(intraday,marketclose)}
  intraday = intraday[2:length(intraday)];
  for (d in 1:cDays) {
    pdatad = pdata[as.character(dates[d])]
    pdatad = aggregatePrice(pdatad, on = on, k = k , marketopen = marketopen, marketclose = marketclose)
    z = xts( rep(1,length(intraday)) , order.by = timeDate( paste(dates[d],as.character(intraday),sep="") , format = "%Y-%m-%d %H:%M:%S"))
    pdatad = merge.xts( z , pdatad )$pdatad
    pdatad = na.locf(pdatad)
    rdatad = makeReturns(pdatad)
    rdatad = rdatad[time(rdatad) > min(time(rdatad))]
    rdata = rbind(rdata, rdatad)
    mR = rbind(mR, as.numeric(rdatad))
  }
  mR[is.na(mR)]=0
  M = ncol(mR)
  if (cDays == 1) {
    mR = as.numeric(rdata)
    estimdailyvol = switch(dailyvol, bipower = rBPCov(mR), 
                           medrv = medRV(mR), rv = RV(mR))
  }else {
    estimdailyvol = switch(dailyvol, bipower = apply(mR, 
                                                     1, "rBPCov"), medrv = apply(mR, 1, "medRV"), rv = apply(mR, 
                                                                                                             1, "RV"))
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
    estimperiodicvol_temp = diurnal(stddata = mstdR, method = periodicvol, 
                                    dummies = dummies, P1 = P1, P2 = P2)[[1]]
    estimperiodicvol = rep(1,M)
    estimperiodicvol[selection] = estimperiodicvol_temp
    mfilteredR = mR/matrix(rep(estimperiodicvol, cDays), 
                           byrow = T, nrow = cDays)
    estimdailyvol = switch(dailyvol, bipower = apply(mfilteredR, 
                                                     1, "rBPCov"), medrv = apply(mfilteredR, 1, "medRV"), 
                           rv = apply(mfilteredR, 1, "RV"))
  }
  out = cbind(rdata, rep(sqrt(estimdailyvol * (1/M)), each = M) * 
                rep(estimperiodicvol, cDays), rep(sqrt(estimdailyvol * 
                                                         (1/M)), each = M), rep(estimperiodicvol, cDays))
  out = xts(out, order.by = time(rdata))
  names(out) = c("returns", "vol", "dailyvol", "periodicvol")
  return(out)
}

#' Stochastic periodicity model
#' 
#' This function estimates the spot volatility by using the stochastic periodcity model of Beltratti & Morana (2001)
#' 
#' @export
stochper =  function(rdata, P=5) 
{
  require(FKF)
  N = ncol(rdata)
  days = nrow(rdata)
  rvector = as.vector(t(rdata)) 
  lambda = (2*pi)/N;
  
  # initial values (let user set them?)
  sigma <- 0.03
  sigma_mu <- 0.005
  sigma_h <- 0.007
  sigma_k <- 0.06
  phi <- 0.194
  rho <- 0.986
  mu <- c(1.87, -0.42)
  delta_c <- c(0.25, -0.05, -0.2, 0.13, 0.02)
  delta_s <- c(-1.2, 0.11, 0.26, -0.03, 0.08)
  
  # transform parameters to allow for unrestricted optimization (domain -Inf to Inf)
  par_t <- c(sigma = log(sigma), sigma_mu = log(sigma_mu), sigma_h = log(sigma_h), sigma_k = log(sigma_k), phi = log(phi/(1-phi)), rho = log(rho/(1-rho)), mu = mu, delta_c = delta_c, delta_s = delta_s)

  out <- optim(par_t, loglikBM, yt = rvector, N = N, days = days, P = P, method="BFGS")
  
  # recreate model to obtain volatility estimates
  ss <- ssmodel(out$par, days, N)
  kf <- fkf(a0 = ss$a0, P0 = ss$P0, dt = ss$dt, ct = ss$ct, Tt = ss$Tt, Zt = ss$Zt, HHt = ss$HHt, GGt = ss$GGt, yt = matrix(rvector, ncol = length(rvector)))
  sigmahat <- exp((ss$Zt%*%kf$at[,1:(N*days)] + ss$ct + 1.27)/2)
  #write.csv(sigmahat, file="sigmahat.csv")
  
  
  estimates <- c(exp(out$par["sigma"]), exp(out$par["sigma_mu"]), exp(out$par["sigma_h"]), exp(out$par["sigma_k"]), exp(out$par["phi"])/(1+exp(out$par["phi"])), exp(out$par["rho"])/(1+exp(out$par["rho"])), out$par[7:18])
  print(estimates)
  return(out)
}

#' Calculate log likelihood using Kalman Filter
#' 
#' This function returns the average log likehood value of the stochastic periodicity model, given the input parameters.
#' 
#'@export
loglikBM <- function(par_t, yt, days, N = 288, P = 5)
{
  ss <- ssmodel(par_t, days, N, P)
  yt = matrix(yt, ncol = length(yt))
  kf <- fkf(a0 = ss$a0, P0 = ss$P0, dt = ss$dt, ct = ss$ct, Tt = ss$Tt, Zt = ss$Zt, HHt = ss$HHt, GGt = ss$GGt, yt = yt)
  #print(kf$logLik/length(yt))
  #print(par_t)
  return(-kf$logLik/length(yt))
}

#' Generate state space model
#' 
#' This function creates the state space matrices from the input parameters.
#' The output is in the format used by the FKF package.
#'
#' @export
ssmodel <- function(par_t, days, N = 288, P = 5)
{
  par <- c(exp(par_t["sigma"]), exp(par_t["sigma_mu"]), exp(par_t["sigma_h"]), exp(par_t["sigma_k"]), exp(par_t["phi"])/(1+exp(par_t["phi"])), exp(par_t["rho"])/(1+exp(par_t["rho"])), par_t[7:18])
  lambda=(2*pi)/288
  a0 <- c(0, 0, par["delta_c1"], par["delta_s1"])
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
  for (k in 2:P)
  {
      c2 <- c2 + par[paste("delta_c", k, sep="")]*cos(k*lambda*n) + par[paste("delta_s", k, sep="")]*sin(k*lambda*n)
  }
  ct <- matrix(ct + c2, ncol = N*days)
  
  return(list(a0 = a0, P0 = P0, Tt = Tt, Zt = Zt, GGt = GGt, HHt = HHt, dt = dt, ct = ct))
}