% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/garch.R
\name{predict.garch1c1}
\alias{predict.garch1c1}
\title{Prediction for GARCH(1,1) time series}
\usage{
\method{predict}{garch1c1}(object, n.ahead = 1, Nsim = 1000, eps,
  sigmasq, seed = NULL, ...)
}
\arguments{
\item{object}{an object from class \code{"garch1c1"}.}

\item{n.ahead}{maximum horizon (lead time) for prediction.}

\item{Nsim}{number of Monte Carlo simulations for simulation based quantities.}

\item{eps}{the time series to predict, only the last value is used.}

\item{sigmasq}{the (squared) volatilities, only the last value is used.}

\item{seed}{an integer, seed for the random number generator.}

\item{...}{currently not used.}
}
\value{
an object from S3 class \code{"predict_garch1c1"} containing
    the following components:
\item{eps}{point predictions (conditional expectations) of the time series (equal
    to zero for pure GARCH).}
\item{h}{point predictions (conditional expectations)of the squared volatilities.}
\item{model}{the model.}
\item{call}{the call.}
\item{pi_plugin}{Prediction intervals for the time series, based on plug-in
    distributions, see Details.}
\item{pi_sim}{Simulation based prediction intervals for the time series, see Details.}
\item{dist_sim}{simulation samples from the predictive distributions of the time
    series and the volatilties. }
}
\description{
Predict GARCH(1,1) time series.
}
\details{
Plug-in prediction intervals and predictive distributions are obtained by
inserting the predicted volatility in the conditional densities. For
predictions more than one lag ahead these are not the real predictive
distributions but the prediction intervals are usually adequate.

For simulation prediction intervals we generate a (large) number of
continuations of the given time series. Prediction intervals can be based on
sample quantiles. The generated samples are stored in the returned object and
can be used for further exploration of the predictive
distributions. \code{dist_sim$eps} contains the simulated future values of
the time series and \code{dist_sim$h} the corresponding (squared)
volatilities.  Both are matrices whose \code{i}-th rows contain the predicted
quantities for horizon \code{i}.
}
\note{
This function is under development and may be changed.
}
\examples{
## set up a model and simulate a time series
mo <- GarchModel(omega = 0.4, alpha = 0.3, beta = 0.5)
a1 <- sim_garch1c1(mo, n = 1000, n.start = 100)

## predictions for T+1,...,T+5 (T = time of last value)
## Nsim is small to reduce the load on CRAN, usually Nsim is larger.
a.pred <- predict(mo, n.ahead = 5, Nsim = 1000, eps = a1$eps, sigmasq = a1$h, seed = 1234)

## preditions for the time series
a.pred$eps

## PI's for eps - plug-in and simulated
a.pred$pi_plugin
a.pred$pi_sim

## a DIY alculation of PI's using the simulated sample paths
t(apply(a.pred$dist_sim$eps, 1, function(x) quantile(x, c(0.025, 0.975))))

## further investigate the predictive distributions
t(apply(a.pred$dist_sim$eps, 1, function(x) summary(x)))

## compare predictive densities for h=2 and h=5
plot(density(a.pred$dist_sim$eps[2, ]), ylim = c(0,.25))
lines(density(a.pred$dist_sim$eps[5, ]), col = "blue")

## predictions of sigma_t^2
a.pred$h

## plug-in predictions of sigma_t
sqrt(a.pred$h)

## simulation predictive densities of sigma_t for h = 2 and h = 5
plot(density(sqrt(a.pred$dist_sim$h[2, ])), xlim = c(0, 6))
lines(density(sqrt(a.pred$dist_sim$h[5, ])), col = "blue")

## VaR and ES for different horizons
cbind(h = 1:5,
      VaR = apply(a.pred$dist_sim$eps, 1, function(x) VaR(x, c(0.05))),
      ES = apply(a.pred$dist_sim$eps, 1, function(x) ES(x, c(0.05))) )

## fit a GARCH(1,1) model to exchange rate data and predict
gmo1 <- fGarch::garchFit(formula = ~garch(1, 1), data = fGarch::dem2gbp,
  include.mean = FALSE, cond.dist = "norm", trace = FALSE)
mocoef <- gmo1@fit$par
mofitted <- GarchModel(omega = mocoef["omega"], alpha = mocoef["alpha1"],
  beta = mocoef["beta1"])
gmo1.pred <- predict(mofitted, n.ahead = 5, Nsim = 1000, eps = gmo1@data,
  sigmasq = gmo1@h.t, seed = 1234)
gmo1.pred$pi_plugin
gmo1.pred$pi_sim

}
