## ----message=FALSE-------------------------------------------------------
library(spotvolatility)
out <- spotvol(sample_prices_5min)
class(out)
sigma_hat <- out$spot

## ------------------------------------------------------------------------
str(sample_prices_5min)
out <- spotvol(sample_prices_5min)

## ------------------------------------------------------------------------
str(sample_returns_5min)
out <- spotvol(sample_returns_5min)

