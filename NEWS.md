# cvar 0.4.1

* changed the JSS reference to use the new-style doi.

* fixed a bug in the tests, in v0.4-0, that was causing faiure of the tests on
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
