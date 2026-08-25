#' Compute Value-at-Risk (VaR)
#'
#' @description
#'
#'     \code{VaR} is a generic function for computation of Value-at-Risk. The default method
#'     computes VaR corresponding to a distribution, usually fitted or implied by a fitted
#'     model.  There are also methods for data (numeric or matrix).
#'
#' @note We use the traditional definition of VaR as the negated lower quantile. For example,
#'     if \eqn{X} are returns on an asset, VAR\eqn{{}_\alpha}{_a} = \eqn{-q_\alpha}{-q_a},
#'     where \eqn{q_\alpha}{-q_a} is the lower \eqn{\alpha}{a} quantile of \eqn{X}.
#'     Equivalently, VAR\eqn{{}_\alpha}{_a} is equal to the lower \eqn{1-\alpha}{1-a}
#'     quantile of \eqn{-X}.
#'
#' @inheritParams ES
#'
#' @param intercept,slope compute VaR for the linear transformation \code{intercept +
#'     slope*X}, where \code{X} has distribution specified by \code{dist}, see Details.
#' 
#'
#' @details
#'
#'     \code{VaR} is S3 generic. The meaning of the parameters for the methods defined in
#'     \pkg{cvar} and described here are the same as in \code{\link{ES}}, including the
#'     recycling rules.
#'
#'     \code{VaR_qf} and \code{VaR_cdf} are streamlined, non-generic, variants for the common
#'     case when the \code{"..."} parameters are scalar. The parameters \code{p_loss},
#'     \code{intercept}, and \code{slope} can be vectors, as for \code{VaR}.
#'
#'     Argument \code{dist} can also be a numeric vector. In that case the ES is computed,
#'     effectively, for the empirical cumulative distribution function (ecdf) of the
#'     vector. The ecdf is not created explicitly and the \code{\link[stats]{quantile}}
#'     function is used instead for the computation of VaR. Arguments in \code{"..."} are
#'     passed eventually to \code{quantile()} and can be used, for example, to select a
#'     non-defult method for the computation of quantiles.
#'
#'     If \code{dist} is of class \code{"matrix"}, then the numeric method is applied to each
#'     column.
#' 
#'     In practice, we may need to compute VaR associated with data. The distribution comes
#'     from fitting a model. In the simplest case, we fit a distribution to the data,
#'     assuming that the sample is i.i.d. For example, a normal distribution \eqn{N(\mu,
#'     \sigma^2)} can be fitted using the sample mean and sample variance as estimates of the
#'     unknown parameters \eqn{\mu} and \eqn{\sigma^2}, see section \sQuote{Examples}. For
#'     other common distributions there are specialised functions to fit their parameters and
#'     if not, general optimisation routines can be used. More soffisticated models may be
#'     used, even time series models such as GARCH and mixture autoregressive models.
#'     
#'     When argument \code{dist} represents log-returns, VaR is for the log-returns. Use
#'     \code{transf = TRUE} to return its value for the returns.
#'
#'     
#' @param tol tollerance
#'
#' @seealso \code{\link{ES}} for ES,
#'
#'     \code{\link[=predict.garch1c1]{predict}} for examples with fitted models
#' 
#' @concept Value at Risk
#'
#' @examples
#' cvar::VaR(qnorm, c(0.01, 0.05), dist.type = "qf")
#'
#' ## the following examples use these values, obtained by fitting a normal distribution to
#' ## some data:
#' muA <- 0.006408553
#' sigma2A <- 0.0004018977
#'
#' ## with quantile function, giving the parameters directly in the call:
#' res1 <- cvar::VaR(qnorm, 0.05, mean = muA, sd = sqrt(sigma2A))
#' res2 <- cvar::VaR(qnorm, 0.05, intercept = muA, slope = sqrt(sigma2A))
#' abs((res2 - res1)) # 0, intercept/slope equivalent to mean/sd
#'
#' ## with quantile function, which already knows the parameters:
#' my_qnorm <- function(p) qnorm(p, mean = muA, sd = sqrt(sigma2A))
#' res1_alt <- cvar::VaR(my_qnorm, 0.05)
#' abs((res1_alt - res1))
#' 
#' ## with cdf the precision depends on solving an equation
#' res1a <- cvar::VaR(pnorm, 0.05, dist.type = "cdf", mean = muA, sd = sqrt(sigma2A))
#' res2a <- cvar::VaR(pnorm, 0.05, dist.type = "cdf", intercept = muA, slope = sqrt(sigma2A))
#' abs((res1a - res2)) # 3.287939e-09
#' abs((res2a - res2)) # 5.331195e-11, intercept/slope better numerically
#'
#' ## as above, but increase the precision, this is probably excessive
#' res1b <- cvar::VaR(pnorm, 0.05, dist.type = "cdf",
#'                    mean = muA, sd = sqrt(sigma2A), tol = .Machine$double.eps^0.75)
#' res2b <- cvar::VaR(pnorm, 0.05, dist.type = "cdf",
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
#'     vAz1 <- cvar::VaR(qnorm, 0.05, mean = mu, sd = sqrt(sigma2))
#'     vAz2 <- cvar::VaR(qnorm, 0.05, intercept = mu, slope = sqrt(sigma2))
#'
#'     vAz1a <- cvar::VaR(pnorm, 0.05, dist.type = "cdf",
#'                        mean = mu, sd = sqrt(sigma2))
#'     vAz2a <- cvar::VaR(pnorm, 0.05, dist.type = "cdf",
#'                        intercept = mu, slope = sqrt(sigma2))
#'
#'     vAz1b <- cvar::VaR(pnorm, 0.05, dist.type = "cdf",
#'                    mean = mu, sd = sqrt(sigma2),
#'                    tol = .Machine$double.eps^0.75)
#'     vAz2b <- cvar::VaR(pnorm, 0.05, dist.type = "cdf",
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
VaR <- function(dist, p_loss = 0.05, ...){
    UseMethod("VaR")
}

warn_about_x <- function(s){
    warning("argument 'x' has been renamed to 'p_loss'. 'x' still works with this version of 'cvar' but will stop doing so in the next release.\n Please change your scripts suitably.\n It most cases there is no need to name this argument when calling ", s)
}

#' @name VaR
#'
#' @export
VaR_qf <- function(dist, p_loss = 0.05, ...,
                   intercept = 0, slope  = 1, tol = .Machine$double.eps^0.5, x){
    if(!missing(x)){
        p_loss <- x
        warn_about_x("VaR_qf")
    }
    dist <- match.fun(dist)
    ## assumes "..." are scalar
    res <-  if(length(p_loss) > 1)
                - sapply(p_loss, dist, ...)
            else
                - dist(p_loss, ...)

    - intercept + slope * res
}

#' @name VaR
#'
#' @export
VaR_cdf <- function(dist, p_loss = 0.05, ...,
                    intercept = 0, slope  = 1, tol = .Machine$double.eps^0.5, x){
    if(!missing(x)){
        p_loss <- x
        warn_about_x("VaR_cdf")
    }
    
    dist <- match.fun(dist)
    ## assumes "..." are scalar
    res <-  if(length(p_loss) > 1)
                - sapply(p_loss, cdf2quantile, ..., MoreArgs = list(cdf = dist, tol = tol))
            else
                - cdf2quantile(p_loss, dist, tol = tol, ...) # TODO: interval?

    - intercept + slope * res
}

#' @name VaR
#'
#' @export
VaR.default <- function(dist, p_loss = 0.05, dist.type = "qf", ...,
                        intercept = 0, slope  = 1, tol = .Machine$double.eps^0.5,
                        transf = FALSE, x){
    if(!missing(x)){
        p_loss <- x
        warn_about_x("VaR")
    }
    
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
            # sapply(p_loss, fu)
    res <- switch(dist.type,
                  qf = - mapply(dist, p_loss, ...),
                  cdf = - mapply(cdf2quantile, p_loss, ..., MoreArgs = list(cdf = dist, tol = tol)),
                  pdf = {
                      stop("Not ready yet. Please supply quantile function or cdf")
                  },
                  stop('argument dist.type should be one of "qf", "cdf" or "pdf"')
                  )
    res <- - intercept + slope * res

    if(transf)
        - (exp(-res) - 1)
    else
        res
}

#' @name VaR
#'
#' @export
VaR.numeric <- function(dist, p_loss = 0.05, ..., intercept = 0, slope  = 1, transf = FALSE, x){
    if(!missing(x)){
        p_loss <- x
        warn_about_x("VaR")
    }
    
    ## dist is raw data here
    res <- - quantile(dist, p_loss, ...)

    res <- - intercept + slope * res
    if(transf)
        - (exp(-res) - 1)
    else
        res
}

#' @name VaR
#'
#' @export
VaR.matrix <- function(dist, p_loss = 0.05, ..., intercept = 0, slope  = 1, transf = FALSE, x){
    if(!missing(x)){
        p_loss <- x
        warn_about_x("VaR")
    }

    lst <- apply(dist, 2, identity, simplify = FALSE) # to use in 'mapply'

    res <- mapply(VaR.numeric, lst, intercept = intercept, slope  = slope, transf = transf,
                  MoreArgs = list(p_loss = p_loss, ...))
    res
}



#' @title Compute expected shortfall (ES)
#'
#' @description
#'     Compute expected shortfall (ES).
#'
#' @details
#' 
#'     \code{ES} is a generic function for computation of expected shortfall. The default
#'     method computes ES corresponding to a distribution, usually fitted or implied by a
#'     fitted model.  There are also methods for data (numeric or matrix).
#'
#'     Alternative terms for ES include CVaR (conditional value at risk), AVaR (average value
#'     at risk), and ETL (expected tail loss).
#' 
#'     The default method of \code{ES} computes the expected shortfall for distributions
#'     specified by the arguments. \code{dist} is typically a function (or the name of
#'     one). What \code{dist} computes is determined by \code{dist.type}, whose default
#'     setting is \code{"qf"} (the quantile function). Other possible settings of
#'     \code{dist.type} include \code{"cdf"} and \code{"pdf"}.  Additional arguments for
#'     \code{dist} can be given with the \code{"..."} arguments.
#'
#'     Argument \code{dist} can also be a numeric vector. In that case the ES is computed,
#'     effectively, for the empirical cumulative distribution function (ecdf) of the
#'     vector. The ecdf is not created explicitly and the \code{\link[stats]{quantile}}
#'     function is used instead for the computation of VaR. Arguments in \code{"..."} are
#'     passed eventually to \code{quantile()} and can be used, for example, to select a
#'     non-defult method for the computation of quantiles.
#'
#'     If \code{dist} is a matrix, the numeric method is applied to each of its columns.
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
#'     A bunch of expected shortfalls is computed if argument \code{p_loss} or any of the
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
#'     When argument \code{dist} represents log-returns, ES is for the log-returns. Use
#'     \code{transf = TRUE} to return its value for the returns. Note that the ES of the
#'     returns cannot be obtained by exponentiating the ES for the log-returns.
#'
#' @param dist specifies the distribution whose ES is computed, usually a function or a name
#'     of a function computing quantiles, cdf, pdf, or a random number generator. Can also be
#'     a numeric vector or matrix, representing data, see section \sQuote{Details}.
#'
#' @param p_loss level, default is 0.05.
#'
#' @param dist.type a character string specifying what is computed by \code{dist}, such as
#'     "qf" or "cdf".
#'
#' @param qf quantile function, only used if \code{dist.type = "pdf"}.
#'
#' @param intercept,slope compute ES for the linear transformation \code{intercept +
#'     slope*X}, where \code{X} has distribution specified by \code{dist}, see Details.
#'
#' @param control additional control parameters for the numerical integration routine.
#'
#' @param transf use only if \code{"dist"} represents log-returns. If \code{TRUE}, give the
#'     result for the returns (not log-returns), see section \sQuote{Details}.
#' 
#' @param ... passed on to \code{dist}.
#'
#' @param x deprecated and will soon be removed. \code{x} was renamed to \code{p_loss},
#'     please use the latter.
#'
#'
#' @return a numeric vector or matrix
#'
#' @seealso \code{\link{VaR}} for VaR,
#'
#'     \code{\link[=predict.garch1c1]{predict}} for examples with fitted models
#' 
#' @concept CVaR
#' @concept AVaR
#' @concept ETL
#' @concept expected shortfall
#' @concept conditional value at risk
#' @concept average value at risk
#' @concept expected tail loss
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
#' ES(pnorm, 0.95, dist.type = "cdf")
#' ES(qnorm, 0.95, dist.type = "qf")
#' ## - VaRES::esnormal(0.95, 0, 1)
#' ## - PerformanceAnalytics::ETL(p=0.05, method = "gaussian", mu = 0,
#' ##                             sigma = 1, weights = 1)             # same
#'
#' cvar::ES(pnorm, dist.type = "cdf")
#' cvar::ES(qnorm, dist.type = "qf")
#' cvar::ES(pnorm, 0.05, dist.type = "cdf")
#' cvar::ES(qnorm, 0.05, dist.type = "qf")
#'
#' ## this uses "pdf"
#' cvar::ES(dnorm, 0.05, dist.type = "pdf", qf = qnorm)
#'
#'
#' ## this gives warning (it does more than simply computing ES):
#' ## PerformanceAnalytics::ETL(p=0.95, method = "gaussian", mu = 0, sigma = 1, weights = 1)
#'
#' ## run this if VaRRES is present
#' \dontrun{
#' x <- seq(0.01, 0.99, length = 100)
#' y <- sapply(x, function(p) cvar::ES(qnorm, p, dist.type = "qf"))
#' yS <- sapply(x, function(p) - VaRES::esnormal(p))
#' plot(x, y)
#' lines(x, yS, col = "blue")
#' }
#'
#' @export
ES <- function(dist, p_loss = 0.05, ...){
    UseMethod("ES")
}

CVaR <- AVaR <- ES # synonyms for ES, currently not exported


#' @name ES
#'
#' @export
ES.default <- function(dist, p_loss = 0.05, dist.type = "qf", qf, ...,
                       intercept = 0, slope = 1, control = list(), transf = FALSE, x){
    if(!missing(x)){
        p_loss <- x
        warn_about_x("ES")
    }
    
    ## TODO: this method needs consolidation; it became repetitive when features were added.
    
    dist <- match.fun(dist)

    stopifnot(length(transf) == 1) # only scalar transf for this method
    if(transf){
        ## can't apply intercept and slope at the end here;
        ## so pass them in each call to VaR or (recursively) ES
        if(length(p_loss) == 1 && all(sapply(list(..., intercept, slope), length) < 2)){
            if(dist.type == "pdf"){
                fu <- function(p) - (exp(-intercept + slope * p) - 1) * dist(p, ...) / p_loss
        
                iargs <- list(f = fu, lower = -Inf, upper = qf(p_loss)) # stop.on.error = FALSE,?
        
            }else{
                ## TODO: check the theoretical justification here
                fu <- function(p){
                    v <- VaR(dist = dist, p_loss = p, dist.type = dist.type, intercept = intercept, slope = slope, ...)
                    - (exp(-v) - 1)  / p_loss
                }
        
                iargs <- list(f = fu, lower = 0, upper = p_loss) # stop.on.error = FALSE,?
            }
        
            iargs[names(control)] <- control
        
            res <- do.call("integrate", iargs)$value
        }else{
            if(dist.type == "pdf"){
                if(is.list(qf))
                    res <- mapply(ES, MoreArgs = list(dist = dist, dist.type = dist.type, control = control, transf = transf, intercept = intercept, slope = slope),
                                  p_loss = p_loss, qf = qf, ...)
                else
                    res <- mapply(ES, MoreArgs = list(dist = dist, dist.type = dist.type, control = control, transf = transf, intercept = intercept, slope = slope, qf = qf),
                              p_loss = p_loss, ...)
                    
            }else{
                res <- mapply(ES, MoreArgs = list(dist = dist, dist.type = dist.type, control = control), intercept = intercept, slope = slope, transf = transf,
                              p_loss = p_loss, ...)
            }
        }
        
        res
    }else{
        if(length(p_loss) == 1 && all(sapply(list(...), length) < 2)){
            if(dist.type == "pdf"){
                fu <- function(p) - p * dist(p, ...) / p_loss
        
                iargs <- list(f = fu, lower = -Inf, upper = qf(p_loss)) # stop.on.error = FALSE,?
        
            }else{
                fu <- function(p) VaR(dist = dist, p_loss = p, dist.type = dist.type, ...) / p_loss
        
                iargs <- list(f = fu, lower = 0, upper = p_loss) # stop.on.error = FALSE,?
            }
        
            iargs[names(control)] <- control
        
            res <- do.call("integrate", iargs)$value
        }else{
            if(dist.type == "pdf"){
                if(is.list(qf))
                    res <- mapply(ES, MoreArgs = list(dist = dist, dist.type = dist.type, control = control),
                                  p_loss = p_loss, qf = qf, ...)
                else
                    res <- mapply(ES, MoreArgs = list(dist = dist, dist.type = dist.type, control = control, qf = qf),
                              p_loss = p_loss, ...)
                    
            }else{
                res <- mapply(ES, MoreArgs = list(dist = dist, dist.type = dist.type, control = control),
                              p_loss = p_loss, ...)
            }
        }
        
        -intercept + slope * res
    }
}

#' @name ES
#'
#' @export
ES.numeric <- function(dist, p_loss = 0.05, ..., intercept = 0, slope = 1, transf = FALSE, x){
    if(!missing(x)){
        p_loss <- x
        warn_about_x("ES")
    }

    v <- VaR.numeric(dist, p_loss = p_loss, ..., intercept = intercept, slope = slope)
    f <- if(transf){
             function(x){
                 logret <- dist[dist <= -x]
                 - mean(exp(logret) - 1)
             }
         }else                          # direct code doesn't work if length(v) > 1):
                                        #     bad <- dist[dist <= - v]
                                        #     res <- - mean(bad)
             function(x) - mean(dist[dist <= -x])
#browser()             
    res <- sapply(v, f)

    res
}

#' @name ES
#'
#' @export
ES.matrix <- function(dist, p_loss = 0.05, ..., intercept = 0, slope = 1, transf = FALSE, x){
    if(!missing(x)){
        p_loss <- x
        warn_about_x("ES")
    }
    
    ## v <- VaR.matrix(dist, p_loss = p_loss, ..., intercept = intercept, slope = slope)

    lst <- apply(dist, 2, identity, simplify = FALSE) # to use in 'mapply'

    res <- mapply(ES.numeric, lst, intercept = intercept, slope  = slope, transf = transf,
                  MoreArgs = list(p_loss = p_loss, ...))
    res
}
