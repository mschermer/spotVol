#' spotvol
#'
#' Estimating Spot Volatility
#' 
#' 
spotvol =  function(pdata, dailyvol = "bipower", periodicvol = "TML", on = "minutes", 
                    k = 5, dummies = FALSE, P1 = 4, P2 = 2,  marketopen = "09:30:00", 
                    marketclose = "16:00:00") 
{
  require(chron);
  require(timeDate);
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