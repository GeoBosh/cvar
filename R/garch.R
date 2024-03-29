#' Specify a GARCH model
#'
#' Specify a GARCH model.
#'
#' Argument \code{model} can be the result of a previous call to \code{GarchModel}.
#' Arguments in \code{"..."} overwrite current components of \code{model}.
#'
#' \code{GarchModel} guarantees that code using it will continue to work
#' transparently for the user even if the internal represedtation of GARCH
#' models in this package is changed or additional functionality is added.
#'
#' @param model a GARCH model or a list.
#' @param ... named arguments specifying the GARCH model.
#' @param model.class a class for the result. By default \code{GarchModel()}
#'     decides the class of the result.
#' @return an object from suitable GARCH-type class
#'
#' @examples
#' ## GARCH(1,1) with Gaussian innovations
#' mo1a <- GarchModel(omega = 1, alpha = 0.3, beta = 0.5)
#' mo1b <- GarchModel(omega = 1, alpha = 0.3, beta = 0.5, cond.dist = "norm")
#'
#' ## equivalently, the parameters can be given as a list
#' p1 <- list(omega = 1, alpha = 0.3, beta = 0.5)
#' mo1a_alt <- GarchModel(p1)
#' mo1b_alt <- GarchModel(p1, cond.dist = "norm")
#' stopifnot(identical(mo1a, mo1a_alt), identical(mo1b, mo1b_alt))
#'
#' ## additional arguments modify values already in 'model'
#' mo_alt <- GarchModel(p1, beta = 0.4)
#'
#' ## set also initial values
#' mo2 <- GarchModel(omega = 1, alpha = 0.3, beta = 0.5, esp0 = - 1.5, h0 = 4.96)
#'
#' ##  GARCH(1,1) with standardised-t_5
#' mot <- GarchModel(omega = 1, alpha = 0.3, beta = 0.5, cond.dist = list("std", nu = 5))
#'
#' @export
GarchModel <- function(model = list(), ..., model.class = NULL){
    ## TODO: check the correctness of the parameters
    ## alpha, beta, cond.dist, 2nd order stationary, initial value for eps_t^2, h_t

    ## 2019-03-15 TODO: handle fGARCH models
    ## if(inherits(model, "fGARCH")){
    ##
    ## }else{

        dots <- list(...)
        if(length(dots) > 0)
            model[names(dots)] <- dots

        class(model) <- "GarchModel0"
        if(is.null(model.class) && length(model$alpha) == 1 && length(model$beta) == 1)
            class(model) <- "garch1c1"                     # GARCH(1,1), 'c' stands for 'comma'
    ## }

    model
}

## 2022-11-01 change .dist from 'list' to 'env', was:
    # .dist <- list(
    #     norm = list(d = call("dnorm", x = NA, mean = 0, sd = 1),
    #                 p = call("pnorm", q = NA, mean = 0, sd = 1),
    #                 q = call("qnorm", p = NA, mean = 0, sd = 1),
    #                 r = call("rnorm", n = NA, mean = 0, sd = 1)),
    #     std  = list(d = call("dstd",  x = NA, mean = 0, sd = 1, nu = NA),
    #                 p = call("pstd",  q = NA, mean = 0, sd = 1, nu = NA),
    #                 q = call("qstd",  p = NA, mean = 0, sd = 1, nu = NA),
    #                 r = call("rstd",  n = NA, mean = 0, sd = 1, nu = NA)),
    #     ged  = list(d = call("dged",  x = NA, mean = 0, sd = 1, nu = NA),
    #                 p = call("pged",  q = NA, mean = 0, sd = 1, nu = NA),
    #                 q = call("qged",  p = NA, mean = 0, sd = 1, nu = NA),
    #                 r = call("rged",  n = NA, mean = 0, sd = 1, nu = NA))
    # )

.dist <- new.env()

## N(0,1) is built-in
.dist$norm <- list(d = call("dnorm", x = NA, mean = 0, sd = 1),
                   p = call("pnorm", q = NA, mean = 0, sd = 1),
                   q = call("qnorm", p = NA, mean = 0, sd = 1),
                   r = call("rnorm", n = NA, mean = 0, sd = 1))


    # generalise and prepare to make fGarch suggested; was:
    # 
    # .dist$std  <- list(d = call("dstd",  x = NA, mean = 0, sd = 1, nu = NA),
    #                    p = call("pstd",  q = NA, mean = 0, sd = 1, nu = NA),
    #                    q = call("qstd",  p = NA, mean = 0, sd = 1, nu = NA),
    #                    r = call("rstd",  n = NA, mean = 0, sd = 1, nu = NA))
    # 
    # .dist$ged  <- list(d = call("dged",  x = NA, mean = 0, sd = 1, nu = NA),
    #                    p = call("pged",  q = NA, mean = 0, sd = 1, nu = NA),
    #                    q = call("qged",  p = NA, mean = 0, sd = 1, nu = NA),
    #                    r = call("rged",  n = NA, mean = 0, sd = 1, nu = NA))

.dist_list <- list(
    ## TODO: allow different packages here by letting 'package be vector or completely
    ##       refactor all this. Only '.get_dist_elem()' should know the implementation.
    std = list(package = "fGarch",
               d = list(str2lang("fGarch::dstd"), x = NA, mean = 0, sd = 1, nu = NA),
               p = list(str2lang("fGarch::pstd"), q = NA, mean = 0, sd = 1, nu = NA),
               q = list(str2lang("fGarch::qstd"), p = NA, mean = 0, sd = 1, nu = NA),
               r = list(str2lang("fGarch::rstd"), n = NA, mean = 0, sd = 1, nu = NA)),

    ged = list(package = "fGarch",
               d = list(str2lang("fGarch::dged"), x = NA, mean = 0, sd = 1, nu = NA),
               p = list(str2lang("fGarch::pged"), q = NA, mean = 0, sd = 1, nu = NA),
               q = list(str2lang("fGarch::qged"), p = NA, mean = 0, sd = 1, nu = NA),
               r = list(str2lang("fGarch::rged"), n = NA, mean = 0, sd = 1, nu = NA))
)

.get_dist_elem <- function(dist, what){
    if(is.list(dist)){
        params <- dist[-1]
        dist <- dist[[1]]
    }else
        params <- list()

    if(is.null(dist)){  ## TODO: this probably should be handled by the caller
        ## normal distribution is default
        res <- .dist$norm[[what]]
    }else{
        stopifnot(is.character(dist), length(dist) == 1)

        gen <- .dist[[dist]]
        if(!is.null(gen))  # already prepared
            res <- gen[[what]]
        else{ # not prepared yet
            gen_list <- .dist_list[[dist]]
            if(is.null(gen_list))
                stop("unknown distribution specified.")

            if(!requireNamespace(gen_list$package))
                stop("requested distribution property needs package ", gen_list$package,
                     "\nplease install package '", gen_list$package, "' and try again")
    
            ## prepare
            ## .dist is env, so assignment is persistent
            .dist[[dist]] <- lapply(gen_list[-1], as.call) # dropping $package
            res <- .dist[[dist]][[what]]
            if(is.null(res))
                stop("property '", "'what not available for '", dist, "'")
        }
    }

    if(length(params) > 0)  # set parameters if specified
        res[names(params)] <- params

    res
}

#' Simulate GARCH(1,1) time series
#'
#' Simulate GARCH(1,1) time series.
#'
#' The simulated time series is in component \code{eps} of the returned value.
#' For exploration of algorithms and eestimation procedures, the volatilities
#' and the standardised innovations are also returned.
#'
#' The random seed at the start of the simulations is saved in the returned
#' object.  A speficific seed can be requested with argument \code{seed}. In
#' that case the simulations are done with the specified seed and the old state
#' of the random number generator is restored before the function returns.
#'
#' @param model a GARCH(1,1) model, an object obtained from \code{GarchModel}.
#' @param n the length of the generated time series.
#' @param n.start number of warm-up values, which are then dropped.
#' @param seed an integer to use for setting the random number generator.
#'
#' @return a list with components:
#' \item{eps}{the time series,}
#' \item{h}{the (squared) volatilities,}
#' \item{eta}{the standardised innovations,}
#' \item{model}{the GARCH(1,1) model,}
#' \item{.sim}{a list containing the parameters of the simulation,}
#' \item{call}{the call.}
#'
#' @note This function is under development and may be changed.
#'
#' @export
sim_garch1c1 <- function(model, n, n.start = 0, seed = NULL){
    ## TODO: seed is not used currently

    garch1c1_fields <- c("omega", "alpha", "beta", "cond.dist", "eps0", "h0", "eps0sq")

    stopifnot(all(names(model) %in% garch1c1_fields))

    omega  <- model$omega
    alpha  <- model$alpha
    beta   <- model$beta
    rgen   <- .get_dist_elem(model$cond.dist, "r")

    eps0   <- model$eps0
    h0     <- model$h0
    eps0sq <- model$eps0sq

    ## see Francq & Zakoian, p. 142, for alternative initialisations
    if(is.null(eps0sq)){
        model$eps0sq <- eps0sq <-
            if(is.null(eps0))
                omega / (1 - alpha - beta)
            else
                eps0^2
    }

    if(is.null(h0)){
        model$h0 <- h0 <-
            omega / (1 - alpha - beta)
    }

    N <- n + n.start
    h <- eps <- numeric(N)

    RNGstate <- .RNGstate(seed)
    ## 2020-03-04 was:    if(!is.null(RNGstate$oldRNGstate)) # or !is.null(seed)
    ##    (see comment at a similar command)
    if(!is.null(seed))
        on.exit(
            if(is.null(RNGstate$RNGstate)){
                ## TRUE id seed is NULL but also if .Random.seed was not set.
                "nothing to do"
            }else if(is.null(RNGstate$oldRNGstate)){
                rm(".Random.seed", envir = .GlobalEnv)
            }else{
                assign(".Random.seed", RNGstate$oldRNGstate, envir = .GlobalEnv)
            }
        )

    rgen$n <- N # or: rgen[[2]] <- N
    eta <- eval(rgen)

    h[1] <- omega + alpha * eps0sq + beta * h0
    for(i in 1:(N - 1)){
        eps[i] <- sqrt(h[i]) * eta[i]
        h[i + 1] <- omega + alpha * eps[i]^2 + beta * h[i]
    }
    eps[N] <- sqrt(h[N]) * eta[N]

    ## TODO: return eta only conditionally?
    ## TODO: include the model in the returned value?
    ## TODO: return also the initialisation values (in a separate component)?
    list(eps = eps[(n.start + 1):N], h = h[(n.start + 1):N], eta = eta[(n.start + 1):N],
         model = model,
         .sim = list(seed = seed, rand.gen = rgen, n = n, n.start = n.start),
         call = match.call()
         )
}

## modelled after the beginning of `stats:::simulate.lm()`
.RNGstate <- function(seed = NULL){
    ## 2020-03-04 - changing the logic below
    ##
    ## if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    ##     runif(1)

    oldRNGstate <- if(exists(".Random.seed", envir = .GlobalEnv))
                       get(".Random.seed", envir = .GlobalEnv)
                   else
                       NULL
    
    if (is.null(seed))
        list(RNGstate = oldRNGstate)
    else{
        ## This doesn't resolve the problem
        ## suppressWarnings(RNGversion("3.5.0"))  # 2019-03-13 temporary, RNG changed in R-devel.
        ##                                        #            see email from Kurt Hornik in Org/
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        ## on.exit(assign(".Random.seed", oldRNGstate, envir = .GlobalEnv))
        list(RNGstate = RNGstate, oldRNGstate = oldRNGstate)
    }
}

#' Prediction for GARCH(1,1) time series
#'
#' Predict GARCH(1,1) time series.
#'
#' Plug-in prediction intervals and predictive distributions are obtained by
#' inserting the predicted volatility in the conditional densities. For
#' predictions more than one lag ahead these are not the real predictive
#' distributions but the prediction intervals are usually adequate.
#'
#' For simulation prediction intervals we generate a (large) number of
#' continuations of the given time series. Prediction intervals can be based on
#' sample quantiles. The generated samples are stored in the returned object and
#' can be used for further exploration of the predictive
#' distributions. \code{dist_sim$eps} contains the simulated future values of
#' the time series and \code{dist_sim$h} the corresponding (squared)
#' volatilities.  Both are matrices whose \code{i}-th rows contain the predicted
#' quantities for horizon \code{i}.
#'
#' The random seed at the start of the simulations is saved in the returned
#' object.  A speficific seed can be requested with argument \code{seed}. In
#' that case the simulations are done with the specified seed and the old state
#' of the random number generator is restored before the function returns.
#' This setup is similar to \code{\link{sim_garch1c1}}.
#'
#' @param object an object from class \code{"garch1c1"}.
#' @param n.ahead maximum horizon (lead time) for prediction.
#' @param Nsim number of Monte Carlo simulations for simulation based quantities.
#' @param eps the time series to predict, only the last value is used.
#' @param sigmasq the (squared) volatilities, only the last value is used.
#' @param seed an integer, seed for the random number generator.
#' @param ... currently not used.
#'
#' @return an object from S3 class \code{"predict_garch1c1"} containing
#'     the following components:
#' \item{eps}{point predictions (conditional expectations) of the time series (equal
#'     to zero for pure GARCH).}
#' \item{h}{point predictions (conditional expectations)of the squared volatilities.}
#' \item{model}{the model.}
#' \item{call}{the call.}
#' \item{pi_plugin}{Prediction intervals for the time series, based on plug-in
#'     distributions, see Details.}
#' \item{pi_sim}{Simulation based prediction intervals for the time series, see Details.}
#' \item{dist_sim}{simulation samples from the predictive distributions of the time
#'     series and the volatilties. }
#'
#' @examples
#' op <- options(digits = 4)
#' 
#' ## set up a model and simulate a time series
#' mo <- GarchModel(omega = 0.4, alpha = 0.3, beta = 0.5)
#' a1 <- sim_garch1c1(mo, n = 1000, n.start = 100, seed = 20220305)
#'
#' ## predictions for T+1,...,T+5 (T = time of last value)
#' ## Nsim is small to reduce the load on CRAN, usually Nsim is larger.
#' a.pred <- predict(mo, n.ahead = 5, Nsim = 1000, eps = a1$eps,
#'   sigmasq = a1$h, seed = 1234)
#'
#' ## preditions for the time series
#' a.pred$eps
#'
#' ## PI's for eps - plug-in and simulated
#' a.pred$pi_plugin
#' a.pred$pi_sim
#'
#' ## a DIY calculation of PI's using the simulated sample paths
#' t(apply(a.pred$dist_sim$eps, 1, function(x) quantile(x, c(0.025, 0.975))))
#'
#' ## further investigate the predictive distributions
#' t(apply(a.pred$dist_sim$eps, 1, function(x) summary(x)))
#'
#' ## compare predictive densities for horizons 2 and 5:
#' h2 <- a.pred$dist_sim$eps[2, ]
#' h5 <- a.pred$dist_sim$eps[5, ]
#' main <- "Predictive densities: horizons 2 (blue) and 5 (black)"
#' plot(density(h5), main = main)
#' lines(density(h2), col = "blue")
#'
#' ## predictions of sigma_t^2
#' a.pred$h
#'
#' ## plug-in predictions of sigma_t
#' sqrt(a.pred$h)
#'
#' ## simulation predictive densities (PD's) of sigma_t for horizons 2 and 5:
#' h2 <- sqrt(a.pred$dist_sim$h[2, ])
#' h5 <- sqrt(a.pred$dist_sim$h[5, ])
#' main <- "PD's of sigma_t for horizons 2 (blue) and 5 (black)"
#' plot(density(h2), col = "blue", main = main)
#' lines(density(h5))
#'
#' ## VaR and ES for different horizons
#' cbind(h = 1:5,
#'       VaR = apply(a.pred$dist_sim$eps, 1, function(x) VaR(x, c(0.05))),
#'       ES = apply(a.pred$dist_sim$eps, 1, function(x) ES(x, c(0.05))) )
#'
#' ## fit a GARCH(1,1) model to exchange rate data and predict
#' gmo1 <- fGarch::garchFit(formula = ~garch(1, 1), data = fGarch::dem2gbp,
#'   include.mean = FALSE, cond.dist = "norm", trace = FALSE)
#' mocoef <- gmo1@fit$par
#' mofitted <- GarchModel(omega = mocoef["omega"], alpha = mocoef["alpha1"],
#'   beta = mocoef["beta1"])
#' gmo1.pred <- predict(mofitted, n.ahead = 5, Nsim = 1000, eps = gmo1@data,
#'   sigmasq = gmo1@h.t, seed = 1234)
#' gmo1.pred$pi_plugin
#' gmo1.pred$pi_sim
#'
#' op <- options(op) # restore options(digits)
#'
#' @note This function is under development and may be changed.
#'
#' @export
predict.garch1c1 <- function(object, n.ahead = 1, Nsim = 1000, eps, sigmasq, seed = NULL, ...){
    model <- object
    omega <- model$omega
    alpha <- model$alpha
    beta  <- model$beta

    model$eps0 <- eps0 <- eps[length(eps)]
    model$h0   <- h0   <- sigmasq[length(sigmasq)] # TODO: compute if not present?

    ## TODO: add non-zero mean and possibly ARMA conditional mean
    pred_eps <- numeric(n.ahead)

    ## tsff slides, p. II.7.4/y1718
    pred_h <- numeric(n.ahead)
    pred_h[1] <- omega + alpha * eps0^2 + beta * h0
    if(n.ahead > 1)
        for(i in 2:n.ahead){
            pred_h[i] <- omega + (alpha + beta) * pred_h[i - 1]
        }

    fq <- .get_dist_elem(model$cond.dist, "q")
    fq$sd <- sqrt(pred_h)

    fq[[2]] <- 0.025
    lwr <- eval(fq)
    fq[[2]] <- 0.975
    upr <- eval(fq)

    eps_pi <- cbind(lwr, upr)

    ## need local variants of withr::with_seed() withr::with_preserve_seed()
    ## but such are not available (as of 2018-09-28), so do it DIY:
    ##
    ## Memorise the state in any case but restore the old state only if 'seed' is not NULL.
    RNGstate <- .RNGstate(seed)
    ## 2020-03-04 was:
    ## 
    ##    if(!is.null(RNGstate$oldRNGstate)) # or !is.null(seed)
    ##
    ## But the check is not equivalent to is.null(seed)!
    ## Moreover, if RNGstate$oldRNGstate is NULL, .Random.seed should be removed,
    ##     (so doing it now)
    if(!is.null(seed))
        on.exit(
            if(is.null(RNGstate$RNGstate)){
                ## TRUE id seed is NULL but also if .Random.seed was not set.
                "nothing to do"
            }else if(is.null(RNGstate$oldRNGstate)){
                rm(".Random.seed", envir = .GlobalEnv)
            }else{
                assign(".Random.seed", RNGstate$oldRNGstate, envir = .GlobalEnv)
            }
        )


    h_sim <- eps_sim <- matrix(NA_real_, nrow = n.ahead, ncol = Nsim)
    for(i in 1:Nsim){
        wrk <- sim_garch1c1(model = model, n = n.ahead, n.start = 0)
        eps_sim[ , i] <- wrk$eps
        h_sim[ , i] <- wrk$h
    }


    sim_pi <- t(apply(eps_sim, 1, function(x) quantile(x, c(0.025, 0.975))))

    structure(list(
        eps = pred_eps, h = pred_h, model = model, call = match.call(),
        pi_plugin = eps_pi, pi_sim = sim_pi,
        dist_sim = list(eps = eps_sim, h = h_sim, RNGstate = RNGstate$RNGstate)
    ), class = "predict_garch1c1")
}
