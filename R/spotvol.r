#' spotvol
#'
#' Estimating Spot Volatility
#' 
#' @export
#' @examples
#' data(sample)
#' vol1 <- spotvol(sample)
#' init = list(sigma = 0.03, sigma_mu = 0.005, sigma_h = 0.007,
#'    sigma_k = 0.06, phi = 0.194, rho = 0.986, mu = c(1.87,-0.42),
#'    delta_c = c(0.25, -0.05, -0.2, 0.13, 0.02), delta_s = c(-1.2, 0.11, 0.26, -0.03, 0.08))
#' vol2 <- spotvol(sample, method = "stochper", init = init)
#' plot(vol1$spotvol, type="l")
#' lines(vol2$spotvol, col="red")
spotvol =  function(pdata, dailyvol = "bipower", periodicvol = "TML", on = "minutes", 
                    k = 5, dummies = FALSE, P1 = 4, P2 = 2,  marketopen = "09:30:00", 
                    marketclose = "16:00:00", method = "detper", init=NULL) 
{
  Sys.setenv(TZ="GMT") # only for convenience during testing
  require(chron);
  require(timeDate);
  require(highfrequency);
  dates = unique(format(time(pdata), "%Y-%m-%d"))
  cDays = length(dates)
  rdata = mR = c()
  if(on=="minutes"){
    intraday = seq(from=times(marketopen), to=times(marketclose), by=times(paste("00:0",k,":00",sep=""))) 
  }
  if(as.character(tail(intraday,1))!=marketclose){intraday=c(intraday,marketclose)}
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

  out = switch(method, 
           detper = detper(mR, cDays, dailyvol = dailyvol, periodicvol = periodicvol, dummies = dummies, P1 = P1, P2 = P2), 
           stochper = stochper(mR, P = max(P1,P2), init = init))
  
  return(out)
}

#' Deterministic periodicity model
#' 
#' From the original spotVol function in highfrequency
detper = function(mR, cDays, dailyvol = "bipower", periodicvol = "TML", dummies = FALSE, P1 = 4, P2 = 2)
{
  M = ncol(mR)
  if (cDays == 1) {
    mR = as.numeric(rdata)
    estimdailyvol = switch(dailyvol, 
                           bipower = rBPCov(mR), 
                           medrv = medRV(mR), rv = RV(mR))
  }else {
    estimdailyvol = switch(dailyvol, 
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
    estimperiodicvol_temp = diurnal(stddata = mstdR, method = periodicvol, 
                                  dummies = dummies, P1 = P1, P2 = P2)[[1]]
    estimperiodicvol = rep(1,M)
    estimperiodicvol[selection] = estimperiodicvol_temp
    mfilteredR = mR/matrix(rep(estimperiodicvol, cDays), 
                         byrow = T, nrow = cDays)
    estimdailyvol = switch(dailyvol, bipower = apply(mfilteredR, 1, "rBPCov"),
                            medrv = apply(mfilteredR, 1, "medRV"), 
                            rv = apply(mfilteredR, 1, "RV"))
    out <- list(spotvol = rep(sqrt(estimdailyvol * (1/M)), each = M) * rep(estimperiodicvol, cDays),
                dailyvol = estimdailyvol,
                periodicvol = estimperiodicvol)
    class(out) <- "spotvol"
    return(out)
  }
}

#' Stochastic periodicity model
#' 
#' This function estimates the spot volatility by using the stochastic periodcity model of Beltratti & Morana (2001)
stochper =  function(mR, init=NULL, P=5) 
{
  require(FKF)
  N = ncol(mR)
  days = nrow(mR)
  
  logr2 = log(mR^2)
  rvector = as.vector(t(logr2)) 
  lambda = (2*pi)/N;
  
  if(is.null(init))
  {
    sigma <- 0.03
    sigma_mu <- 0.005
    sigma_h <- 0.007
    sigma_k <- 0.06
    phi <- 0.194
    rho <- 0.986
    mu <- c(1.87, -0.42)
    delta_c <- c(0.25, -0.05, -0.2, 0.13, 0.02)
    delta_s <- c(-1.2, 0.11, 0.26, -0.03, 0.08) 
    
    par_t <- c(sigma = log(sigma), sigma_mu = log(sigma_mu), sigma_h = log(sigma_h), sigma_k = log(sigma_k),
               phi = log(phi/(1-phi)), rho = log(rho/(1-rho)), mu = mu, delta_c = delta_c, delta_s = delta_s)
  }
  else
  {
    par_t <- c(sigma = log(init$sigma), sigma_mu = log(init$sigma_mu), sigma_h = log(init$sigma_h),
               sigma_k = log(init$sigma_k), phi = log(init$phi/(1-init$phi)), rho = log(init$rho/(1-init$rho)),
               mu = init$mu, delta_c = init$delta_c, delta_s = init$delta_s)   
  }
  
  # transform parameters to allow for unrestricted optimization (domain -Inf to Inf)

  opt <- optim(par_t, loglikBM, yt = rvector, N = N, days = days, P = P, method="BFGS")
  
  # recreate model to obtain volatility estimates
  ss <- ssmodel(opt$par, days, N)
  kf <- fkf(a0 = ss$a0, P0 = ss$P0, dt = ss$dt, ct = ss$ct, Tt = ss$Tt, Zt = ss$Zt, 
            HHt = ss$HHt, GGt = ss$GGt, yt = matrix(rvector, ncol = length(rvector)))
  sigmahat <- as.vector(exp((ss$Zt%*%kf$at[,1:(N*days)] + ss$ct + 1.27)/2))
  
  estimates <- c(exp(opt$par["sigma"]), exp(opt$par["sigma_mu"]), exp(opt$par["sigma_h"]), exp(opt$par["sigma_k"]),
                 exp(opt$par["phi"])/(1+exp(opt$par["phi"])), exp(opt$par["rho"])/(1+exp(opt$par["rho"])), opt$par[7:18])

  out <- list(spotvol = sigmahat,
              par = estimates)
  class(out) <- "spotvol"
  return(out)
}

#' Calculate log likelihood using Kalman Filter
#' 
#' This function returns the average log likehood value of the stochastic periodicity model, given the input parameters.
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
ssmodel <- function(par_t, days, N = 288, P = 5)
{
  par <- c(exp(par_t["sigma"]), exp(par_t["sigma_mu"]), exp(par_t["sigma_h"]), exp(par_t["sigma_k"]), 
           exp(par_t["phi"])/(1+exp(par_t["phi"])), exp(par_t["rho"])/(1+exp(par_t["rho"])), par_t[7:18])
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
      plot(seas, main = "Non-parametric and parametric periodicity estimates", 
           xlab = "intraday period", type = "l", lty = 3)
      legend("topright", c("Parametric", "Non-parametric"), cex = 1.1,
             lty = c(1,3), lwd = 1, bty = "n")
      seas = highfrequency:::diurnalfit(theta = theta, P1 = P1, P2 = P2, intraT = intraT, 
                                        dummies = dummies)
      lines(seas, lty = 1)
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

#' Sample data
#' 
#' 'sample_real5minprices' from highfrequency package
#' 
#' @docType data
#' @name sample
#' @usage data(sample)
#' @format xts
NULL