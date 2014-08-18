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

## ----results='hide'---------------------------------------
out1 <- spotvol(sample_prices_5min, method = "detper")
out2 <- spotvol(sample_prices_5min, method = "stochper")

## ---------------------------------------------------------
names(out1)
names(out2)

## ---------------------------------------------------------
vol_detper <- spotvol(sample_returns_5min, method = "detper")
plot(vol_detper)

## ----crop=TRUE,fig.height=4-------------------------------
plot(as.numeric(t(vol_detper$spot[1:3,])), type = "l")

## ----eval=FALSE-------------------------------------------
#  plot(vol_detper, length = 3*288)

## ----tidy=TRUE--------------------------------------------
init <- list(sigma = 0.03, sigma_mu = 0.005, sigma_h = 0.007, sigma_k = 0.06, phi = 0.194, rho = 0.986, mu = c(1.87,-0.42), delta_c = c(0.25, -0.05, -0.2, 0.13, 0.02), delta_s = c(-1.2, 0.11, 0.26, -0.03, 0.08))
vol_stochper <- spotvol(sample_returns_5min, method = "stochper", init = init, control = list(maxit = 20))

## ----crop=TRUE,fig.height=4-------------------------------
plot(as.numeric(t(vol_detper$spot[1:3,])), type = "l")
lines(as.numeric(t(vol_stochper$spot[1:3,])), col = "red")

## ---------------------------------------------------------
h1 = bw.nrd0((1:nrow(sample_returns_5min))*(5*60))
vol3 <- spotvol(sample_returns_5min, method = "kernel", h = h1)

## ----crop=TRUE,fig.height=3,tidy=TRUE---------------------
quarticity = (nrow(sample_returns_5min)/3)*rowSums(sample_returns_5min^4)
vol4 <- spotvol(sample_returns_5min, method = "kernel", est = "quarticity")
plot(log(quarticity), vol4$par$h)

## ----results='hide'---------------------------------------
vol5 <- spotvol(sample_returns_5min, method = "kernel", est = "cv")

## ----crop=TRUE,fig.height=4,tidy=TRUE---------------------
plot(vol3, length = 2880)
lines(as.numeric(t(vol4$spot))[1:2880], col="red")
lines(as.numeric(t(vol5$spot))[1:2880], col="blue")
legend("topright", c("h = simple estimate", "h = quarticity corrected",
  "h = crossvalidated"), col = c("black", "red", "blue"), lty=1)

## ----tidy=TRUE--------------------------------------------
simdata <- matrix(sqrt(5/3)*rt(3000, df = 5), ncol = 500, byrow = TRUE)
simdata <- c(1, 1, 1.5, 1.5, 2, 1)*simdata

## ----crop=TRUE,fig.height=3,tidy=TRUE---------------------
vol6 <- spotvol(simdata, method = "piecewise", m = 200, n = 100, online = FALSE)
plot(vol6)

## ----crop=TRUE,fig.height=3,tidy=TRUE---------------------
vol7 <- spotvol(simdata, method = "piecewise", m = 200, n = 100, online = TRUE, volest = "tau")
plot(vol7)

## ----crop=TRUE,fig.height=3,tidy=TRUE---------------------
vol8 <- spotvol(sample_returns_5min, method = "garch", model = "sGARCH", solver.control = list(maxeval = 1000))
vol9 <- spotvol(sample_returns_5min, method = "garch", model = "eGARCH", solver.control = list(maxeval = 1000))
plot(as.numeric(t(vol8$spot))[6000:7000], type = "l")
lines(as.numeric(t(vol9$spot))[6000:7000], col = "red")
legend("topleft", c("GARCH", "eGARCH"), col = c("black", "red"), lty=1)

