spotVol
=======

Repository for R package spotvolatility, which was developed during the Google Summer of Code 2014

Overview: https://www.google-melange.com/gsoc/project/details/google/gsoc2014/maarten473/5668600916475904

To install the package, run the following R commands:

library(devtools)

install_github("mschermer/spotVol")

Description:

The spotvolatility package offers several methods to estimate spot volatility and its intraday seasonality, using high-frequency data.

The following spot volatility estimation methods have been implemented:

Deterministic periodicity

The spot volatility is decomposed into a a deterministic periodic factor f_i (identical for every day in the sample) and a daily factor s_t (identical for all observations within a day). Both components are then estimated separately. For more details, see Taylor and Xu (1997) and Andersen and Bollerslev (1997). The jump robust versions by Boudt et al. (2011) have also been implemented.

Stochastic periodicity

This method by Beltratti and Morana (2001) assumes the periodicity factor to be stochastic. The spot volatility estimation is split into four components: a random walk, an autoregressive process, a stochastic cyclical process and a deterministic cyclical process. The model is estimated using a quasi-maximum likelihood method based on the Kalman Filter. The package FKF is used to apply the Kalman filter.

Nonparametric filtering

This method by Kristensen (2010) filters the spot volatility in a nonparametric way by applying kernel weights to the standard realized volatility estimator. Different kernels and bandwidths can be used to focus on specific characteristics of the volatility process.

Piecewise constant volatility

Another nonparametric method is that of Fried (2012), which assumes the volatility to be piecewise constant over local windows. Robust two-sample tests are applied to detect changes in variability between subsequent windows. The spot volatility can then be estimated by evaluating regular realized volatility estimators within each local window.

GARCH models with intraday seasonality

The package also includes an option to apply GARCH models, implemented by the rugarch package, to estimate spot volatility from intraday data. This is done by including external regressors in the model. These regressors are based on a flexible Fourier form, which was also used in the stochastic and deterministic periodicity estimation methods.

References

Andersen, T. G. and T. Bollerslev (1997). Intraday periodicity and volatility persistence in financial markets. Journal of Empirical Finance 4, 115-158.

Beltratti, A. and C. Morana (2001). Deterministic and stochastic methods for estimation of intraday seasonal components with high frequency data. Economic Notes 30, 205-234.

Boudt K., Croux C. and Laurent S. (2011). Robust estimation of intraweek periodicity in volatility and jump detection. Journal of Empirical Finance 18, 353-367.

Fried, Roland (2012). On the online estimation of local constant volatilities. Computational Statistics and Data Analysis 56, 3080-3090.

Kristensen, Dennis (2010). Nonparametric filtering of the realized spot volatility: A kernel-based approach. Econometric Theory 26, 60-93.

Taylor, S. J. and X. Xu (1997). The incremental volatility information in one million foreign exchange quotations. Journal of Empirical Finance 4, 317-340.


