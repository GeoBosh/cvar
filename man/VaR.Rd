% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VaR.R
\name{VaR}
\alias{VaR}
\alias{VaR_qf}
\alias{VaR_cdf}
\alias{VaR.default}
\alias{VaR.numeric}
\title{Compute Value-at-Risk (VaR)}
\usage{
VaR(dist, x = 0.05, dist.type = "qf", ..., intercept = 0,
  slope = 1, tol = .Machine$double.eps^0.5)

VaR_qf(dist, x = 0.05, ..., intercept = 0, slope = 1,
  tol = .Machine$double.eps^0.5)

VaR_cdf(dist, x = 0.05, ..., intercept = 0, slope = 1,
  tol = .Machine$double.eps^0.5)

\method{VaR}{default}(dist, x = 0.05, dist.type = "qf", ...,
  intercept = 0, slope = 1, tol = .Machine$double.eps^0.5)

\method{VaR}{numeric}(dist, x = 0.05, ..., intercept = 0, slope = 1)
}
\arguments{
\item{dist}{specifies the distribution whose ES is computed, usually a function or a name
of a function computing quantiles, cdf, pdf, or a random number generator, see
Details.}

\item{x}{level, default is 0.05}

\item{dist.type}{a character string specifying what is computed by \code{dist}, such as
"qf" or "cdf".}

\item{...}{passed on to \code{dist}.}

\item{intercept}{requests the ES for the linear transformation \code{intercept +
slope*X}, where \code{X} has distribution specified by \code{dist}, see Details.}

\item{slope}{requests the ES for the linear transformation \code{intercept +
slope*X}, where \code{X} has distribution specified by \code{dist}, see Details.}

\item{tol}{tollerance}
}
\description{
\code{Var} computes the Value-at-Risk of the distribution specified by the
    arguments. The meaning of the parameters is the same as in \code{\link{ES}}, including
    the recycling rules.
}
\details{
\code{VaR} is S3 generic. The meaning of the parameters for its default method is the
    same as in \code{\link{ES}}, including the recycling rules.

    \code{VaR_qf} and \code{VaR_cdf} are streamlined, non-generic, variants for the common
    case when the \code{"..."} parameters are scalar. The parameters \code{x},
    \code{intercept}, and \code{slope} can be vectors, as for \code{VaR}.
}
\note{
We use the traditional definition of VaR as the negated lower quantile. For example,
    if \eqn{X} are returns on an asset, VAR\eqn{{}_\alpha}{_a} = \eqn{-q_\alpha}{-q_a},
    where \eqn{q_\alpha}{-q_a} is the lower \eqn{\alpha}{a} quantile of \eqn{X}.
    Equivalently, VAR\eqn{{}_\alpha}{_a} is equal to the lower \eqn{1-\alpha}{1-a}
    quantile of \eqn{-X}.
}
\examples{
cvar::VaR(qnorm, x = c(0.01, 0.05), dist.type = "qf")

## the following examples use these values:
muA <- 0.006408553
sigma2A <- 0.0004018977

## with quantile function
res1 <- cvar::VaR(qnorm, x = 0.05, mean = muA, sd = sqrt(sigma2A))
res2 <- cvar::VaR(qnorm, x = 0.05, intercept = muA, slope = sqrt(sigma2A))
abs((res2 - res1)) # 0, intercept/slope equivalent to mean/sd

## with cdf the precision depends on solving an equation
res1a <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf", mean = muA, sd = sqrt(sigma2A))
res2a <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf", intercept = muA, slope = sqrt(sigma2A))
abs((res1a - res2)) # 3.287939e-09
abs((res2a - res2)) # 5.331195e-11, intercept/slope better numerically

## as above, but increase the precision, this is probably excessive
res1b <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf",
                   mean = muA, sd = sqrt(sigma2A), tol = .Machine$double.eps^0.75)
res2b <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf",
                   intercept = muA, slope = sqrt(sigma2A), tol = .Machine$double.eps^0.75)
abs((res1b - res2)) # 6.938894e-18 # both within machine precision
abs((res2b - res2)) # 1.040834e-16

## relative precision is also good
abs((res1b - res2)/res2) # 2.6119e-16 # both within machine precision
abs((res2b - res2)/res2) # 3.91785e-15


## an extended example with vector args, if "PerformanceAnalytics" is present
if (requireNamespace("PerformanceAnalytics", quietly = TRUE)) withAutoprint({
    data(edhec, package = "PerformanceAnalytics")
    mu <- apply(edhec, 2, mean)
    sigma2 <- apply(edhec, 2, var)
    musigma2 <- cbind(mu, sigma2)

    ## compute in 2 ways with cvar::VaR
    vAz1 <- cvar::VaR(qnorm, x = 0.05, mean = mu, sd = sqrt(sigma2))
    vAz2 <- cvar::VaR(qnorm, x = 0.05, intercept = mu, slope = sqrt(sigma2))

    vAz1a <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf",
                       mean = mu, sd = sqrt(sigma2))
    vAz2a <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf",
                       intercept = mu, slope = sqrt(sigma2))

    vAz1b <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf",
                   mean = mu, sd = sqrt(sigma2),
                   tol = .Machine$double.eps^0.75)
    vAz2b <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf",
                   intercept = mu, slope = sqrt(sigma2),
                   tol = .Machine$double.eps^0.75)

    ## analogous calc. with PerformanceAnalytics::VaR
    vPA <- apply(musigma2, 1, function(x)
        PerformanceAnalytics::VaR(p = .95, method = "gaussian", invert = FALSE,
                                  mu = x[1], sigma = x[2], weights = 1))
    ## the results are numerically the same
    max(abs((vPA - vAz1))) # 5.551115e-17
    max(abs((vPA - vAz2))) #   ""

    max(abs((vPA - vAz1a))) # 3.287941e-09
    max(abs((vPA - vAz2a))) #  1.465251e-10, intercept/slope better

    max(abs((vPA - vAz1b))) # 4.374869e-13
    max(abs((vPA - vAz2b))) # 3.330669e-16
})

}
