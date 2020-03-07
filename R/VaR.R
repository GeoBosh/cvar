#' Compute Value-at-Risk (VaR)
#'
#' @description \code{Var} computes the Value-at-Risk of the distribution specified by the
#'     arguments. The meaning of the parameters is the same as in \code{\link{ES}}, including
#'     the recycling rules.
#'
#' @note We use the traditional definition of VaR as the negated lower quantile. For example,
#'     if \eqn{X} are returns on an asset, VAR\eqn{{}_\alpha}{_a} = \eqn{-q_\alpha}{-q_a},
#'     where \eqn{q_\alpha}{-q_a} is the lower \eqn{\alpha}{a} quantile of \eqn{X}.
#'     Equivalently, VAR\eqn{{}_\alpha}{_a} is equal to the lower \eqn{1-\alpha}{1-a}
#'     quantile of \eqn{-X}.
#'
#' @inheritParams ES
#'
#' @details
#'     \code{VaR} is S3 generic. The meaning of the parameters for its default method is the
#'     same as in \code{\link{ES}}, including the recycling rules.
#'
#'     \code{VaR_qf} and \code{VaR_cdf} are streamlined, non-generic, variants for the common
#'     case when the \code{"..."} parameters are scalar. The parameters \code{x},
#'     \code{intercept}, and \code{slope} can be vectors, as for \code{VaR}.
#'
#' @param tol tollerance
#'
#' @examples
#' cvar::VaR(qnorm, x = c(0.01, 0.05), dist.type = "qf")
#'
#' ## the following examples use these values:
#' muA <- 0.006408553
#' sigma2A <- 0.0004018977
#'
#' ## with quantile function
#' res1 <- cvar::VaR(qnorm, x = 0.05, mean = muA, sd = sqrt(sigma2A))
#' res2 <- cvar::VaR(qnorm, x = 0.05, intercept = muA, slope = sqrt(sigma2A))
#' abs((res2 - res1)) # 0, intercept/slope equivalent to mean/sd
#'
#' ## with cdf the precision depends on solving an equation
#' res1a <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf", mean = muA, sd = sqrt(sigma2A))
#' res2a <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf", intercept = muA, slope = sqrt(sigma2A))
#' abs((res1a - res2)) # 3.287939e-09
#' abs((res2a - res2)) # 5.331195e-11, intercept/slope better numerically
#'
#' ## as above, but increase the precision, this is probably excessive
#' res1b <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf",
#'                    mean = muA, sd = sqrt(sigma2A), tol = .Machine$double.eps^0.75)
#' res2b <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf",
#'                    intercept = muA, slope = sqrt(sigma2A), tol = .Machine$double.eps^0.75)
#' abs((res1b - res2)) # 6.938894e-18 # both within machine precision
#' abs((res2b - res2)) # 1.040834e-16
#'
#' ## relative precision is also good
#' abs((res1b - res2)/res2) # 2.6119e-16 # both within machine precision
#' abs((res2b - res2)/res2) # 3.91785e-15
#'
#'
#' ## an extended example with vector args, if "PerformanceAnalytics" is present
#' if (requireNamespace("PerformanceAnalytics", quietly = TRUE)) withAutoprint({
#'     data(edhec, package = "PerformanceAnalytics")
#'     mu <- apply(edhec, 2, mean)
#'     sigma2 <- apply(edhec, 2, var)
#'     musigma2 <- cbind(mu, sigma2)
#'
#'     ## compute in 2 ways with cvar::VaR
#'     vAz1 <- cvar::VaR(qnorm, x = 0.05, mean = mu, sd = sqrt(sigma2))
#'     vAz2 <- cvar::VaR(qnorm, x = 0.05, intercept = mu, slope = sqrt(sigma2))
#'
#'     vAz1a <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf",
#'                        mean = mu, sd = sqrt(sigma2))
#'     vAz2a <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf",
#'                        intercept = mu, slope = sqrt(sigma2))
#'
#'     vAz1b <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf",
#'                    mean = mu, sd = sqrt(sigma2),
#'                    tol = .Machine$double.eps^0.75)
#'     vAz2b <- cvar::VaR(pnorm, x = 0.05, dist.type = "cdf",
#'                    intercept = mu, slope = sqrt(sigma2),
#'                    tol = .Machine$double.eps^0.75)
#'
#'     ## analogous calc. with PerformanceAnalytics::VaR
#'     vPA <- apply(musigma2, 1, function(x)
#'         PerformanceAnalytics::VaR(p = .95, method = "gaussian", invert = FALSE,
#'                                   mu = x[1], sigma = x[2], weights = 1))
#'     ## the results are numerically the same
#'     max(abs((vPA - vAz1))) # 5.551115e-17
#'     max(abs((vPA - vAz2))) #   ""
#'
#'     max(abs((vPA - vAz1a))) # 3.287941e-09
#'     max(abs((vPA - vAz2a))) #  1.465251e-10, intercept/slope better
#'
#'     max(abs((vPA - vAz1b))) # 4.374869e-13
#'     max(abs((vPA - vAz2b))) # 3.330669e-16
#' })
#'
#' @export
VaR <- function(dist, x = 0.05, dist.type = "qf", ...,
                intercept = 0, slope  = 1, tol = .Machine$double.eps^0.5){
    UseMethod("VaR")
}

#' @name VaR
#'
#' @export
VaR_qf <- function(dist, x = 0.05, ...,
                   intercept = 0, slope  = 1, tol = .Machine$double.eps^0.5){
    dist <- match.fun(dist)
    ## assumes "..." are scalar
    res <-  if(length(x) > 1)
                - sapply(x, dist, ...)
            else
                - dist(x, ...)

    - intercept + slope * res
}

#' @name VaR
#'
#' @export
VaR_cdf <- function(dist, x = 0.05, ...,
                    intercept = 0, slope  = 1, tol = .Machine$double.eps^0.5){
    dist <- match.fun(dist)
    ## assumes "..." are scalar
    res <-  if(length(x) > 1)
                - sapply(x, cdf2quantile, ..., MoreArgs = list(cdf = dist, tol = tol))
            else
                - cdf2quantile(x, dist, tol = tol, ...) # TODO: interval?

    - intercept + slope * res
}

#' @name VaR
#'
#' @export
VaR.default <- function(dist, x = 0.05, dist.type = "qf", ...,
                        intercept = 0, slope  = 1, tol = .Machine$double.eps^0.5){
    dist <- match.fun(dist)
    ## TODO: add "rand" to dist.type
    ## TODO: include 'dist' in the condition and mapply?
            # fu <- switch(dist.type,
            #              qf = function(y) - dist(y, ...),
            #              cdf = function(y) - cdf2quantile(y, dist, tol = tol, ...),
            #              pdf = {
            #                  stop("Not ready yet. Please supply quantile function or cdf")
            #              },
            #              stop('argument dist.type should be one of "qf", "cdf" or "pdf"')
            #              )
            # sapply(x, fu)
    res <- switch(dist.type,
                  qf = - mapply(dist, x, ...),
                  cdf = - mapply(cdf2quantile, x, ..., MoreArgs = list(cdf = dist, tol = tol)),
                  pdf = {
                      stop("Not ready yet. Please supply quantile function or cdf")
                  },
                  stop('argument dist.type should be one of "qf", "cdf" or "pdf"')
                  )
    - intercept + slope * res
}

#' @name VaR
#'
#' @export
VaR.numeric <- function(dist, x = 0.05, ..., intercept = 0, slope  = 1){
    ## dist is raw data here
    res <- - quantile(dist, x, ...)

    - intercept + slope * res
}



## VaR <- function(dist, x = 0.05, dist.type = "qf", ...,
##                 intercept = 0, slope  = 1, tol = .Machine$double.eps^0.5){
##     dist <- match.fun(dist)
##     ## TODO: add "rand" to dist.type
##     ## TODO: include 'dist' in the condition and mapply()?
##     if(length(x) == 1 && all(sapply(list(...), length) < 2))
##         res <- switch(dist.type,
##                       qf = - dist(x, ...),
##                       cdf = - cdf2quantile(x, dist, tol = tol, ...), # TODO: interval?
##                       pdf = {
##                           stop("Not ready yet. Please supply quantile function or cdf")
##                       },
##                       stop('argument dist.type should be one of "qf", "cdf" or "pdf"')
##                       )
##     else{
##             # fu <- switch(dist.type,
##             #              qf = function(y) - dist(y, ...),
##             #              cdf = function(y) - cdf2quantile(y, dist, tol = tol, ...),
##             #              pdf = {
##             #                  stop("Not ready yet. Please supply quantile function or cdf")
##             #              },
##             #              stop('argument dist.type should be one of "qf", "cdf" or "pdf"')
##             #              )
##             # sapply(x, fu)
##         res <- - switch(dist.type,
##                         qf = mapply(dist, x, ...),
##                         cdf = mapply(cdf2quantile, x, ..., MoreArgs = list(cdf = dist, tol = tol)),
##                         pdf = {
##                             stop("Not ready yet. Please supply quantile function or cdf")
##                         },
##                         stop('argument dist.type should be one of "qf", "cdf" or "pdf"')
##                         )
##
##     }
##     - intercept + slope * res
## }

#' @title Compute expected shortfall (ES) of distributions
#'
#' @description \code{ES} computes the expected shortfall for distributions specified by the
#'     arguments. \code{dist} is typically a function (or the name of one). What \code{dist}
#'     computes is determined by \code{dist.type}, whose default setting is \code{"qf"} (the
#'     quantile function). Other possible settings of \code{dist.type} include \code{"cdf"}
#'     and \code{"pdf"}.  Additional arguments for \code{dist} can be given with the
#'     \code{"..."} arguments.
#'
#'     Except for the exceptions discussed below, a function computing VaR for the specified
#'     distribution is constructed and the expected shortfall is computed by numerically
#'     integrating it. The numerical integration can be fine-tuned with argument
#'     \code{control}, which should be a named list, see \code{\link{integrate}} for the
#'     available options.
#'
#'     If \code{dist.type} is \code{"pdf"}, VaR is not computed, Instead, the partial
#'     expectation of the lower tail is computed by numerical integration of \code{x *
#'     pdf(x)}.  Currently the quantile function is required anyway, via argument \code{qf},
#'     to compute the upper limit of the integral. So, this case is mainly for testing and
#'     comparison purposes.
#'
#'
#'     A bunch of expected shortfalls is computed if argument \code{x} or any of the
#'     arguments in \code{"..."} are of length greater than one. They are recycled to equal
#'     length, if necessary, using the normal \R recycling rules.
#'
#'     \code{intercept} and \code{slope} can be used to compute the expected shortfall for
#'     the location-scale transformation \code{Y = intercept + slope * X}, where the
#'     distribution of \code{X} is as specified by the other parameters and \code{Y} is the
#'     variable of interest. The expected shortfall of \code{X} is calculated and then
#'     transformed to that of \code{Y}. Note that the distribution of \code{X} doesn't need
#'     to be standardised, although it typically will.
#'
#'     The \code{intercept} and the \code{slope} can be vectors. Using them may be
#'     particularly useful for cheap calculations in, for example, forecasting, where the
#'     predictive distributions are often from the same family, but with different location
#'     and scale parameters. Conceptually, the described treatment of \code{intercept} and
#'     \code{slope} is equivalent to recycling them along with the other arguments, but more
#'     efficiently.
#'
#'     The names, \code{intercept} and \code{slope}, for the location and scale parameters
#'     were chosen for their expressiveness and to minimise the possibility for a clash with
#'     parameters of \code{dist} (e.g., the Gamma distribution has parameter \code{scale}).
#'
#'
#' @param dist specifies the distribution whose ES is computed, usually a function or a name
#'     of a function computing quantiles, cdf, pdf, or a random number generator, see
#'     Details.
#'
#' @param x level, default is 0.05
#'
#' @param dist.type a character string specifying what is computed by \code{dist}, such as
#'     "qf" or "cdf".
#'
#' @param qf quantile function, only used if \code{dist.type = "pdf"}.
#'
#' @param intercept,slope requests the ES for the linear transformation \code{intercept +
#'     slope*X}, where \code{X} has distribution specified by \code{dist}, see Details.
#'
#' @param control additional control parameters for the numerical integration routine.
#'
#' @param ... passed on to \code{dist}.
#'
#'
#' @return a numeric vector
#'
#' @examples
#' ES(qnorm)
#'
#' ## Gaussian
#' ES(qnorm, dist.type = "qf")
#' ES(pnorm, dist.type = "cdf")
#'
#' ## t-dist
#' ES(qt, dist.type = "qf", df = 4)
#' ES(pt, dist.type = "cdf", df = 4)
#'
#' ES(pnorm, x= 0.95, dist.type = "cdf")
#' ES(qnorm, x= 0.95, dist.type = "qf")
#' ## - VaRES::esnormal(0.95, 0, 1)
#' ## - PerformanceAnalytics::ETL(p=0.05, method = "gaussian", mu = 0,
#' ##                             sigma = 1, weights = 1)             # same
#'
#' cvar::ES(pnorm, dist.type = "cdf")
#' cvar::ES(qnorm, dist.type = "qf")
#' cvar::ES(pnorm, x= 0.05, dist.type = "cdf")
#' cvar::ES(qnorm, x= 0.05, dist.type = "qf")
#'
#' ## this uses "pdf"
#' cvar::ES(dnorm, x = 0.05, dist.type = "pdf", qf = qnorm)
#'
#'
#' ## this gives warning (it does more than simply computing ES):
#' ## PerformanceAnalytics::ETL(p=0.95, method = "gaussian", mu = 0, sigma = 1, weights = 1)
#'
#' ## run this if VaRRES is present
#' \dontrun{
#' x <- seq(0.01, 0.99, length = 100)
#' y <- sapply(x, function(p) cvar::ES(qnorm, x = p, dist.type = "qf"))
#' yS <- sapply(x, function(p) - VaRES::esnormal(p))
#' plot(x, y)
#' lines(x, yS, col = "blue")
#' }
#'
#' @export
ES <- function(dist, x = 0.05, dist.type = "qf", qf, ...,
               intercept = 0, slope = 1, control = list()){

    ## 2018-09-30 handle "numeric" to complement VaR.numeric()
    if(is.numeric(dist)){
        v <- VaR.numeric(dist, x = 0.05, ..., intercept = intercept, slope = slope)
        bad <- dist[dist <= - v]
        res <- - mean(bad)
        return(res)
    }

    dist <- match.fun(dist)

    if(length(x) == 1 && all(sapply(list(...), length) < 2)){
        if(dist.type == "pdf"){
            fu <- function(p) - p * dist(p, ...) / x

            iargs <- list(f = fu, lower = -Inf, upper = qf(x)) # stop.on.error = FALSE,?

        }else{
            fu <- function(p) VaR(dist = dist, x = p, dist.type = dist.type, ...) / x

            iargs <- list(f = fu, lower = 0, upper = x) # stop.on.error = FALSE,?
        }

        iargs[names(control)] <- control

        res <- do.call("integrate", iargs)$value
    }else{
        ## TODO: this is lazy
        if(dist.type == "pdf"){
            if(is.list(qf))
                res <- mapply(ES, MoreArgs = list(dist = dist, dist.type = dist.type, control = control),
                              x = x, qf = qf, ...)
            else
                res <- mapply(ES, MoreArgs = list(dist = dist, dist.type = dist.type, control = control, qf = qf),
                          x = x, ...)
                
        }else{
            res <- mapply(ES, MoreArgs = list(dist = dist, dist.type = dist.type, control = control),
                          x = x, ...)
        }
    }

    -intercept + slope * res
}

CVaR <- AVaR <- ES
