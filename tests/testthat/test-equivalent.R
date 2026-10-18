set.seed(21)
ultra <- ape::rcoal(25)
ultra$edge.length <- ultra$edge.length / max(ape::node.depth.edgelength(ultra))
fossil <- ape::rphylo(25, birth = 1, death = 0.6, fossils = TRUE)

test_that("OU with a fixed root equals accelerating early burst on ultrametric trees", {
  eq <- kernel_equivalent(k_ou(root = "fixed"), k_eb(upper = 50), ultra,
                          pars = list(ou1 = list(alpha = 1.5)))
  expect_equal(eq$verdict, "identical")
  expect_lt(eq$KL, 1e-6)
  expect_equal(eq[["k2:eb1.rate"]], 3, tolerance = 1e-3)      # r = 2 alpha
})

test_that("a random root is equivalent only when the mean is unknown", {
  p <- list(ou1 = list(alpha = 2))
  free <- kernel_equivalent(k_ou(), k_eb(upper = 50), ultra, pars = p)
  zero <- kernel_equivalent(k_ou(), k_eb(upper = 50), ultra, pars = p, mean = "zero")
  expect_equal(free$verdict, "identical")
  expect_gt(zero$KL, 1e-3)
})

test_that("fossil tips separate OU from early burst", {
  H <- max(ape::node.depth.edgelength(fossil))
  eq <- kernel_equivalent(k_ou(root = "fixed"), k_eb(upper = 100), fossil,
                          pars = list(ou1 = list(alpha = 8 / H)))
  expect_equal(eq$verdict, "distinguishable")
  expect_false(eq$k2_at_bound)
})

test_that("Pagel's lambda equals BM plus white noise on ultrametric trees", {
  eq <- kernel_equivalent(k_lambda(), k_bm() + k_white(), ultra, n_grid = 3)
  expect_true(all(eq$verdict == "identical"))
})

test_that("bounds that prevent a match are flagged", {
  eq <- kernel_equivalent(k_ou(root = "fixed"), k_eb(upper = 1), ultra,
                          pars = list(ou1 = list(alpha = 5)))
  expect_true(eq$k2_at_bound)
  expect_output(print(eq), "at their bounds")
})

test_that("kernels with latent values are rejected", {
  expect_error(kernel_equivalent(k_relaxed(), k_bm(), ultra), "latent")
})
