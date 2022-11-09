# cvar 0.5 (CRAN)

* made `ES` generic (`VaR` was already generic).

* moved `fGarch` from Imports to Suggests.

* renamed argument `x` of `VaR` and `ES` to `p_loss`. `p_loss` seems more
  expressive and suggests that it relates to the losses, usually small numbers
  like `0.05`. Other suitable names like `alpha`, `p`, and `prob`, are commonly
  used as arguments to other functions that might be used as argument `dist` and
  make them more difficult to pass via the `...` arguments.

  For now, an warning is issued if `x` is used as a named argument in a call
  (e.g. `VaR(dist, x = 0.05)`) with the intend to turn that in an error in the
  next release of the package. This change should not be noticed by most users
  since it is much more natural not to name this argument and use something like
  `VaR(dist, 0.05)`.

* moved `fGarch` from Imports to Suggests.


# cvar 0.4.1 (CRAN)

* when the input was numeric, `ES()` was not handling the level `x` properly
  (fixes issue #2, reported by Marius Bommert).

* changed the JSS reference to use the new-style doi.

* fixed a bug in the tests, in v0.4-0, that was causing failure of the tests on
  travis, despite all checks on CRAN passing with OK. `devtools::test()` was
  failing too, but only on the first run in a session, details in the git
  commit.

* set up GHA.


# cvar 0.4-0 (CRAN)

* fix tests to pass with the changed R random generator.

* some new examples and minor documentation changes.


# cvar 0.3-0 (CRAN)

* now `\VignetteIndexEntry` in `Guide_cvar.Rnw` is plain text.

* added experimental support for GARCH models - currently GARCH(1,1) (the API
  may change).

* now the first argument of `VaR()` and `ES()` can be a numeric vector. This is
  useful, e.g., for computing VaR by simulation.

* bugfix:  in `VaR_cdf()` and `VaR_qf()`,  the code for the `if/else` clauses
  had been wrongly swappped. 


# cvar 0.2-0 (CRAN)

* suggest 'covr'.

* setup for travisCI and Coveralls.

* more tests.

* corrected \VignetteIndexEntry in the vignette --- I used the vignette for Rdpack as a
  template but didn't change this entry, which resulted in the vignette appearing with a
  puzzling title on CRAN and other sites.!

* added author@R in DESCRIPTION.

* added experimental web site (docs/)


# cvar 0.1-1 (CRAN)

* added a doi to DESCRIPTION.


# cvar 0.1-0

* first public release.
