## ----setup, include=FALSE---------------------------------
options(width=60)
library(knitr)
knit_hooks$set(crop=function(before, options, envir){
if (before) par(mar=c(2,2,1,1))
})

## ----eval=FALSE-------------------------------------------
#  library(devtools)
#  install_github("mschermer/spotVol")

## ----eval=FALSE-------------------------------------------
#  library(tools)
#  library(spotvolatility)
#  buildVignettes("spotvolatility")

## ----message=FALSE----------------------------------------
library(spotvolatility)
out <- spotvol(sample_prices_5min)
class(out)
sigma_hat <- out$spot

## ----tidy=TRUE--------------------------------------------
str(sample_prices_5min)
out <- spotvol(sample_prices_5min)

## ---------------------------------------------------------
str(sample_returns_5min)
out <- spotvol(sample_returns_5min)

