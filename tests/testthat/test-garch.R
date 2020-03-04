library(cvar)
context("garch")

test_that("GarchModel works ok", {
    mo1a <- GarchModel(omega = 1, alpha = 0.3, beta = 0.5)
    expect_equal(class(mo1a), "garch1c1")

    mo1b <- GarchModel(omega = 1, alpha = 0.3, beta = 0.5, cond.dist = "norm")

    identical(GarchModel(mo1a, omega = 0.4),
              GarchModel(      omega = 0.4, alpha = 0.3, beta = 0.5) )

    ##  GARCH(1,1) with standardised-t_5
    mo_t <- GarchModel(omega = 1, alpha = 0.3, beta = 0.5, cond.dist = list("std", nu = 5))
    expect_identical(mo_t, GarchModel(mo1a, cond.dist = list("std", nu = 5)))
})

test_that("garch1c1 related functions work ok", {
    ## try to deal with the misterious error on TravisCI and devtools::test()
    ## set.seed(123)
    
    a_mo <- GarchModel(omega = 0.4, alpha = 0.3, beta = 0.5)
    a <- sim_garch1c1(a_mo, n = 100, n.start = 100, seed = 1234)
    a_pred <- predict(a_mo, n.ahead = 5, Nsim = 100, eps = a$eps, sigmasq = a$h, seed = 1235)

    ## 2019-03-13 deal with change in RNG in R-devel (for 3.6.0), see email from Kurt Hornik in Org/
    ##        I couldn't make his suggestion work.
    ## TODO: maybe it would be better to just load 'a' from a saved version.
    if(getRversion() < "3.6.0"){
        a_saved <- "a_before_6.0.RDS"
        a_pred_saved <- "a_pred_before_6.0.RDS"
    }else{
        a_saved <- "a.RDS"
        a_pred_saved <- "a_pred.RDS"
    }

    expect_equal_to_reference(a, a_saved)
    expect_equal_to_reference(a_pred, a_pred_saved)

    ## as above but without 'seed'
    sim_garch1c1(a_mo, n = 100, n.start = 100)
    predict(a_mo, n.ahead = 5, Nsim = 100, eps = a$eps, sigmasq = a$h)

    expect_equal(.rgen(NULL), .dist$norm$r)
    ## for now
    expect_error(.rgen("unknowndist"))
    expect_error(.rgen(list("unknowndist")))
    expect_error(.rgen(5))

    .rgen(list("norm"))
    .rgen(list("norm", mean = 5, sd = 3))
    .rgen(list("norm", mean = 5, sd = 3, n = 10))

    .rgen("std")
    .rgen(list("std", df = 5))
    .rgen(list("std", df = 5, n = 10))

    .rgen(list("ged"))

    expect_equal(.get_cond_dist(NULL, "p"), .dist$norm[["p"]])

    .get_cond_dist("norm", "d")
    .get_cond_dist("norm", "p")
    .get_cond_dist("norm", "q")
    .get_cond_dist("norm", "r")

    .get_cond_dist(list("norm"), "p")
    .get_cond_dist(list("norm", mean = 4, sd = 2), "p")
    
    ## for now
    expect_error(.get_cond_dist("unknowndist"))
    expect_error(.get_cond_dist(list("unknowndist")))
    expect_error(.get_cond_dist(5))

})
