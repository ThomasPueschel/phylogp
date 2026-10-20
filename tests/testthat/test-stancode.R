test_that("branch-wise kernels use the pruning likelihood", {
  code <- kernel_stancode(k_lambda() + k_white())
  expect_match(code, "tree_pruning_lpdf(r |", fixed = TRUE)
  expect_match(code, "tree_loo_moments", fixed = TRUE)
  expect_false(grepl("matrix[N, E] B;", code, fixed = TRUE))
})

test_that("dense components switch to dense algebra", {
  code <- kernel_stancode(k_bm() + k_dist(phi = "exponential"))
  expect_match(code, "multi_normal_cholesky_lpdf", fixed = TRUE)
  expect_match(code, "gaussian_loo_lpd", fixed = TRUE)
})

test_that("OU is branch-wise on ultrametric trees and dense otherwise", {
  set.seed(5)
  expect_match(kernel_stancode(k_ou(), ape::rcoal(10)), "tree_pruning_lpdf(r |", fixed = TRUE)
  expect_match(kernel_stancode(k_ou(), ape::rtree(10)), "multi_normal_cholesky_lpdf", fixed = TRUE)
})

test_that("relaxed rates integrate the left-out tip's rate", {
  code <- kernel_stancode(k_relaxed() + k_white())
  expect_match(code, "log_sum_exp(lp)", fixed = TRUE)
  dense <- kernel_stancode(k_relaxed() + k_dist())
  expect_match(dense, "gaussian_loo_moments", fixed = TRUE)
  expect_error(kernel_stancode(k_relaxed() + k_white(df = 4)), "At most one")
})

test_that("standata carries bounds and prior scales for every component", {
  set.seed(6)
  tr <- ape::rcoal(12)
  d <- kernel_standata(k_lambda() + k_white(), tr, y = rnorm(12))
  expect_equal(d$lambda1_lambda_ub, 1)
  expect_true(all(c("lambda1_sigma_scale", "white1_sigma_scale") %in% names(d)))
})
