

# cvar

[![CRANStatusBadge](http://www.r-pkg.org/badges/version/cvar)](https://cran.r-project.org/package=cvar)
[![Build Status](https://travis-ci.com/GeoBosh/cvar.svg?branch=master)](https://travis-ci.com/GeoBosh/cvar)


## Overview

Compute expected shortfall (ES) and Value at Risk (VaR) from a
quantile function, distribution function, random number generator or
probability density function.  ES is also known as Conditional Value
at Risk (CVaR). Virtually any continuous distribution can be
specified.  The functions are vectorised over the arguments.
The computations are done directly from the definitions, see e.g. Acerbi
and Tasche (2002).


## Installing cvar

The [latest stable version](https://cran.r-project.org/package=cvar) is on CRAN. 

    install_packages("cvar")

The vignette shipping with the package gives illustrative examples.
`vignette("Guide_cvar", package = "cvar")`.
(Due to a mixed-up index entry, it appears to have a puzzling title on the [CRAN cvar page](https://cran.r-project.org/package=cvar),
but clicking on it brings up the correct vignette.)

You can install the [development version](https://github.com/GeoBosh/cvar) of `cvar` from Github:

    library(devtools)
    install_github("GeoBosh/cvar")


## Overview

Package `cvar` is a small `R` package with, essentially two
functions &#x2014; `ES` for computing the expected shortfall
and `VaR` for Value at Risk.  The user specifies the
distribution by supplying one of the functions that define a
continuous distribution&#x2014;currently this can be a quantile
function (qf), cumulative distribution function (cdf) or
probability density function (pdf). Virtually any continuous
distribution can be specified.

The functions are vectorised over the parameters of the
distributions, making bulk computations more convenient, for
example for forecasting or model evaluation.

The name of this package, "cvar", comes from *Conditional Value at
Risk* (CVaR), which is an alternative term for expected shortfall.

We chose to use the standard names `ES` and `VaR`,
despite the possibility for name clashes with same named
functions in other packages, rather than invent possibly
difficult to remember alternatives. Just call the functions as
`cvar::ES` and `cvar::VaR` if necessary.

Locations-scale transformations can be specified separately
from the other distribution parameters. This is useful when
such parameters are not provided directly by the distribution
at hand. The use of these parameters often leads to more
efficient computations and better numerical accuracy even if
the distribution has its own parameters for this purpose. Some
of the examples for `VaR` and `ES` illustrate this
for the Gaussian distribution.

Since VaR is a quantile, functions computing it for a given
distribution are convenience functions. `VaR` exported by
`cvar` could be attractive in certain workflows because of
its vectorised distribution parameters, the location-scale
transformation and the possibility to compute it from cdf's
when quantile functions are not available.

