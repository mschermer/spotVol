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

