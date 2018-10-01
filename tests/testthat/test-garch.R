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
    a_mo <- GarchModel(omega = 0.4, alpha = 0.3, beta = 0.5)
    a <- sim_garch1c1(a_mo, n = 100, n.start = 100, seed = 1234)
    a_pred <- predict(a_mo, n.ahead = 5, Nsim = 100, eps = a$eps, sigmasq = a$h, seed = 1235)

    expect_equal_to_reference(a, "a.RDS")
    expect_equal_to_reference(a_pred, "a_pred.RDS")
})
