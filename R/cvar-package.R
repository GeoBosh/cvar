#' @author Georgi N. Boshnakov
#'
#' @importFrom stats quantile runif
#' @importFrom gbutils cdf2quantile
#' @importFrom Rdpack reprompt
#'
#' @importFrom fGarch dstd pstd qstd rstd dged pged qged rged
#'
#' @docType package
#' @name cvar-package
#' @aliases cvar
#'
#' @title Compute Conditional Value-at-Risk and Value-at-Risk
#'
#' @description Compute expected shortfall (ES) and Value at Risk (VaR)
#'     from a quantile function, distribution function, random number
#'     generator or probability density function.  ES is also known as
#'     Conditional Value at Risk (CVaR). Virtually any continuous
#'     distribution can be specified.  The functions are vectorised over
#'     the arguments. Some support for GARCH models is provided, as well.
#'
#' @details
#'
#'     There is a huge number of functions for computations with
#'     distributions in core \R and in contributed packages. Pdf's,
#'     cdf's, quantile functions and random number generators are
#'     covered comprehensively. The coverage of expected shortfall is
#'     more patchy but a large collection of distributions, including
#'     functions for expected shortfall, is provided by
#'     \insertCite{VaRES2013;textual}{cvar}.
#'     \insertCite{PerformanceAnalytics2018;textual}{cvar} and
#'     \insertCite{actuarJSS2008;textual}{cvar} provide packages
#'     covering comprehensively various aspects of risk measurement,
#'     including some functions for expected shortfall.
#'
#'     Package \pkg{cvar} is a small package with, essentially, two main
#'     functions --- \code{ES} for computing the expected shortfall
#'     and \code{VaR} for Value at Risk.  The user specifies the
#'     distribution by supplying one of the functions that define a
#'     continuous distribution---currently this can be a quantile
#'     function (qf), cumulative distribution function (cdf) or
#'     probability density function (pdf). Virtually any continuous
#'     distribution can be specified.
#'
#'     The functions are vectorised over the parameters of the
#'     distributions, making bulk computations more convenient, for
#'     example for forecasting or model evaluation.
#'
#'     The name of this package, "cvar", comes from \emph{Conditional
#'     Value at Risk} (CVaR), which is an alternative term for
#'     expected shortfall.
#'
#'     We chose to use the standard names \code{ES} and \code{VaR},
#'     despite the possibility for name clashes with same named
#'     functions in other packages, rather than invent possibly
#'     difficult to remember alternatives. Just call the functions as
#'     \code{cvar::ES} and \code{cvar::VaR} if necessary.
#'
#'     Locations-scale transformations can be specified separately
#'     from the other distribution parameters. This is useful when
#'     such parameters are not provided directly by the distribution
#'     at hand. The use of these parameters often leads to more
#'     efficient computations and better numerical accuracy even if
#'     the distribution has its own parameters for this purpose. Some
#'     of the examples for \code{VaR} and \code{ES} illustrate this
#'     for the Gaussian distribution.
#'
#'     Since VaR is a quantile, functions computing it for a given
#'     distribution are convenience functions. \code{VaR} exported by
#'     \pkg{cvar} could be attractive in certain workflows because of
#'     its vectorised distribution parameters, the location-scale
#'     transformation, and the possibility to compute it from cdf's
#'     when quantile functions are not available.
#'
#'     Some support for GARCH models is provided, as well. It is
#'     currently under development, see \code{\link{predict.garch1c1}}
#'     for current functionality.
#'
#'     In practice, we may need to compute VaR associated with data. The distribution comes
#'     from fitting a model. In the simplest case, we fit a distribution to the data,
#'     assuming that the sample is i.i.d. For example, a normal distribution \eqn{N(\mu,
#'     \sigma^2)} can be fitted using the sample mean and sample variance as estimates of the
#'     unknown parameters \eqn{\mu} and \eqn{\sigma^2}, see section \sQuote{Examples}. For other
#'     common distributions there are specialised functions to fit their parameters and if
#'     not, general optimisation routines can be used. More soffisticated models may be used,
#'     even time series models such as GARCH and mixture autoregressive models.
#'
#' @references
#'   \insertAllCited{}
#'
#' @seealso
#'   \code{\link{ES}},
#'   \code{\link{VaR}}
#' 
#' @examples
#' ## see the examples for ES(), VaR(), predict.garch1c1()
#'
NULL

.onLoad <- function(lib, pkg){
    Rdpack::Rdpack_bibstyles(package = pkg, authors = "LongNames")
    invisible(NULL)
}
