set.seed(4)
ultra <- ape::rcoal(40)
fossil <- ape::rtree(40)

test_that("any distance kernel is valid on an ultrametric tree", {
  for (phi in c("exponential", "matern32", "matern52", "sqexp")) {
    expect_true(all(kernel_check(k_dist(phi = phi), ultra)$valid))
  }
})

test_that("only completely monotone distance kernels are guaranteed on other trees", {
  expect_true(all(kernel_check(k_dist(phi = "exponential"), fossil)$valid))
  expect_true(all(kernel_check(k_ou(), fossil)$valid))
  chk <- kernel_check(k_dist(phi = "sqexp"), fossil)
  expect_false(chk$valid)
  expect_lt(chk$min_value, 0)
})

test_that("user-defined edge kernels are checked for negative edge variances", {
  ok <- k_edge("pow(len, kappa) .* exp(beta * he)",
               pars = list(kappa = c(0, 1), beta = c(-5, 5)))
  bad <- k_edge("len + beta * he", pars = list(beta = c(-5, 5)))
  expect_true(all(kernel_check(ok, fossil)$valid))
  expect_false(all(kernel_check(bad, fossil)$valid))
})

test_that("kernel_plot returns the correlation matrix", {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  R <- kernel_plot(k_bm() + k_white(), ultra)
  grDevices::dev.off()
  expect_equal(dim(R), c(40, 40))
  expect_equal(unname(diag(R)), rep(1, 40))
})
