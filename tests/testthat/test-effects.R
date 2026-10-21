# plot_effects(): the expected response along one predictor

test_that("plot_effects() returns x'beta on a grid, for numeric and factor predictors", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "CmdStan not installed")
  set.seed(5)
  tr <- ape::rcoal(40)
  x <- kernel_simulate(k_bm(), tr)[, 1]
  g <- factor(sample(c("a", "b"), 40, TRUE))
  d <- data.frame(x = x, g = g, row.names = names(x))
  d$y <- 1 + 0.5 * d$x + 0.3 * (d$g == "b") + kernel_simulate(k_bm() + k_white(), tr)[rownames(d), 1]
  fit <- suppressMessages(phylo_gp(y ~ x + g, d, tr, k_bm() + k_white(), method = "optimize"))
  e <- plot_effects(fit, n = 5, plot = FALSE)
  expect_equal(sort(unique(e$effect)), c("g", "x"))
  b <- coef(fit)
  ex <- e[e$effect == "x", ]
  expect_equal(ex$estimate, b[[1]] + b[[2]] * ex$value + 0, tolerance = 1e-8)  # g at reference level
  eg <- e[e$effect == "g", ]
  expect_equal(eg$level, c("a", "b"))
  expect_equal(eg$estimate, b[[1]] + b[[2]] * mean(d$x) + c(0, b[[3]]), tolerance = 1e-8)
  expect_true(all(is.na(e$lower)))                     # no band for maximum likelihood
  expect_error(plot_effects(fit, "z", plot = FALSE), "Not a predictor")
})
