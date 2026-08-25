library(cvar)
context("VaR")

test_that("VaR works ok", {

muA <- 0.006408553
sigma2A <- 0.0004018977

## with quantile function
res1 <- cvar::VaR(qnorm, p_loss = 0.05, mean = muA, sd = sqrt(sigma2A))
res2 <- cvar::VaR(qnorm, p_loss = 0.05, intercept = muA, slope = sqrt(sigma2A))
abs((res2 - res1)) # 0, intercept/slope equivalent to mean/sd

expect_error(cvar::VaR(dnorm, dist.type = "pdf", p_loss = 0.05, intercept = muA, slope = sqrt(sigma2A)), "Not ready")
cvar::VaR_qf(qnorm, p_loss = 0.05, mean = muA, sd = sqrt(sigma2A))

## with cdf the precision depends on solving an equation
res1a <- cvar::VaR(pnorm, p_loss = 0.05, dist.type = "cdf", mean = muA, sd = sqrt(sigma2A))
res2a <- cvar::VaR(pnorm, p_loss = 0.05, dist.type = "cdf", intercept = muA, slope = sqrt(sigma2A))
abs((res1a - res2)) # 3.287939e-09
abs((res2a - res2)) # 5.331195e-11, intercept/slope better numerically

expect_equal(res1, VaR_qf(qnorm, p_loss = 0.05, mean = muA, sd = sqrt(sigma2A)))
expect_equal(res2a, VaR_cdf(pnorm, p_loss = 0.05, mean = muA, sd = sqrt(sigma2A)))

 
set.seed(1236)
a.num <- rnorm(100)
VaR(a.num)

## test the fix for issue #2
expect_true(ES(a.num) == ES(a.num, p_loss = 0.05))
expect_true(ES(a.num) != ES(a.num, p_loss = 0.01))

## as above, but increase the precision, this is probably excessive
res1b <- cvar::VaR(pnorm, p_loss = 0.05, dist.type = "cdf", mean = muA, sd = sqrt(sigma2A), tol = .Machine$double.eps^0.75)
res2b <- cvar::VaR(pnorm, p_loss = 0.05, dist.type = "cdf", intercept = muA, slope = sqrt(sigma2A), tol = .Machine$double.eps^0.75)

expect_lt( abs((res1b - res2)), 1e-12) # 6.938894e-18 # both within machine precision
expect_lt( abs((res2b - res2)), 1e-12) # 1.040834e-16

## relative precision is also good
abs((res1b - res2)/res2) # 2.6119e-16 # both within machine precision
abs((res2b - res2)/res2) # 3.91785e-15


## an extended example with vector args, if "PerformanceAnalytics" is present
if(require("PerformanceAnalytics")){
    data(edhec, package = "PerformanceAnalytics")
    mu <- apply(edhec, 2, mean)
    sigma2 <- apply(edhec, 2, var)
    musigma2 <- cbind(mu, sigma2)

    ## compute in 2 ways with cvar::VaR
    vAz1 <- cvar::VaR(qnorm, p_loss = 0.05, mean = mu, sd = sqrt(sigma2))
    vAz2 <- cvar::VaR(qnorm, p_loss = 0.05, intercept = mu, slope = sqrt(sigma2))

    vAz1a <- cvar::VaR(pnorm, p_loss = 0.05, dist.type = "cdf",
                       mean = mu, sd = sqrt(sigma2))
    vAz2a <- cvar::VaR(pnorm, p_loss = 0.05, dist.type = "cdf",
                       intercept = mu, slope = sqrt(sigma2))

    vAz1b <- cvar::VaR(pnorm, p_loss = 0.05, dist.type = "cdf",
                   mean = mu, sd = sqrt(sigma2),
                   tol = .Machine$double.eps^0.75)
    vAz2b <- cvar::VaR(pnorm, p_loss = 0.05, dist.type = "cdf",
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

    ## compute in 2 ways with cvar::VaR; this time with transf = TRUE
    ## TODO: need more meaningful example(s) for this
    vAz1tr <- cvar::VaR(qnorm, p_loss = 0.05, mean = mu, sd = sqrt(sigma2), transf = TRUE)
    vAz2tr <- cvar::VaR(qnorm, p_loss = 0.05, intercept = mu, slope = sqrt(sigma2), transf = TRUE)

    vAz1atr <- cvar::VaR(pnorm, p_loss = 0.05, dist.type = "cdf",
                       mean = mu, sd = sqrt(sigma2), transf = TRUE)
    vAz2atr <- cvar::VaR(pnorm, p_loss = 0.05, dist.type = "cdf",
                       intercept = mu, slope = sqrt(sigma2), transf = TRUE)

    vAz1btr <- cvar::VaR(pnorm, p_loss = 0.05, dist.type = "cdf",
                   mean = mu, sd = sqrt(sigma2),
                   tol = .Machine$double.eps^0.75, transf = TRUE)
    vAz2btr <- cvar::VaR(pnorm, p_loss = 0.05, dist.type = "cdf",
                   intercept = mu, slope = sqrt(sigma2),
                   tol = .Machine$double.eps^0.75, transf = TRUE)
}

## the matrix method
m20 <- matrix(rnorm(40), ncol = 2)
colnames(m20) <- c("A", "B")
expect_equal(VaR(m20), VaR(m20, 0.05))
VaR(m20)
VaR(m20, c(0.01, 0.05))

## ES
## use 500 values to mamke 1% and 5% quantiles more likely to be different
m500 <- matrix(rnorm(250), ncol = 2)
ES(m500)
ES(m500, 0.01)
ES(m500, 0.05)
ES(m500, c(0.01, 0.05))

## transf = TRUE
ES(m500, transf = TRUE)
ES(m500, 0.01, transf = TRUE)
ES(m500, 0.05, transf = TRUE)
ES(m500, c(0.01, 0.05), transf = TRUE)

VaR(m500, transf = TRUE)
VaR(m500, 0.01, transf = TRUE)
VaR(m500, 0.05, transf = TRUE)
VaR(m500, c(0.01, 0.05), transf = TRUE)

})

test_that("ES works ok", {

expect_equal(ES(qnorm, dist.type = "qf"),
             ES(qnorm) )

## Gaussian
expect_equal(ES(qnorm, dist.type = "qf"),
             ES(pnorm, dist.type = "cdf") )


## t-dist
expect_equal(ES(qt, dist.type = "qf", df = 4),
             ES(pt, dist.type = "cdf", df = 4) )

expect_equal(ES(pnorm, p_loss = 0.95, dist.type = "cdf"),
             ES(qnorm, p_loss = 0.95, dist.type = "qf") )
## - VaRES::esnormal(0.95, 0, 1)
## - PerformanceAnalytics::ETL(p=0.05, method = "gaussian", mu = 0,
##                             sigma = 1, weights = 1)             # same


## this uses "pdf"
expect_equal(ES(dnorm, p_loss = 0.05, dist.type = "pdf", qf = qnorm),
             ES(pnorm, p_loss = 0.05, dist.type = "cdf") )

## several in one call
ES(pnorm, p_loss = c(0.1, 0.05, 0.01), dist.type = "cdf")
ES(dnorm, p_loss = c(0.1, 0.05, 0.01), dist.type = "pdf", qf = qnorm)
ES(dnorm, p_loss = 0.05, dist.type = "pdf", qf = qnorm)
ES(dnorm, p_loss = c(0.1, 0.05, 0.01), dist.type = "pdf", qf = qnorm)
ES(dnorm, p_loss = c(0.1, 0.05, 0.01), dist.type = "pdf", qf = list(qnorm, qnorm, qnorm))
ES(dnorm, p_loss = c(0.1, 0.05, 0.01), dist.type = "pdf", qf = c(qnorm, qnorm, qnorm))


#### the above with transf = TRUE; TODO: need more meaningful examples
####
expect_equal(ES(qnorm, dist.type = "qf", transf = TRUE),
             ES(qnorm, transf = TRUE) )

## Gaussian
expect_equal(ES(qnorm, dist.type = "qf", transf = TRUE),
             ES(pnorm, dist.type = "cdf", transf = TRUE) )


## t-dist
expect_lt(abs(ES(qt, dist.type = "qf", df = 4, transf = TRUE) - 
                 ES(pt, dist.type = "cdf", df = 4, transf = TRUE)), 1e-7 )

expect_equal(ES(pnorm, p_loss = 0.95, dist.type = "cdf", transf = TRUE),
             ES(qnorm, p_loss = 0.95, dist.type = "qf", transf = TRUE) )
## - VaRES::esnormal(0.95, 0, 1)
## - PerformanceAnalytics::ETL(p=0.05, method = "gaussian", mu = 0,
##                             sigma = 1, weights = 1)             # same


## this uses "pdf"
expect_lt(abs(ES(dnorm, p_loss = 0.05, dist.type = "pdf", qf = qnorm, transf = TRUE) -
              ES(pnorm, p_loss = 0.05, dist.type = "cdf", transf = TRUE)), 1e-7 )

## several in one call
expect_lt(max(abs(
    ES(pnorm, p_loss = c(0.1, 0.05, 0.01), dist.type = "cdf", transf = TRUE) -
    ES(dnorm, p_loss = c(0.1, 0.05, 0.01), dist.type = "pdf", qf = qnorm, transf = TRUE)
    )), 1e-7)

ES(dnorm, p_loss = 0.05, dist.type = "pdf", qf = qnorm, transf = TRUE)
ES(dnorm, p_loss = c(0.1, 0.05, 0.01), dist.type = "pdf", qf = qnorm, transf = TRUE)

expect_equal(
    ES(dnorm, p_loss = c(0.1, 0.05, 0.01), dist.type = "pdf", qf = list(qnorm, qnorm, qnorm), transf = TRUE),
    ES(dnorm, p_loss = c(0.1, 0.05, 0.01), dist.type = "pdf", qf = c(qnorm, qnorm, qnorm), transf = TRUE))



})

