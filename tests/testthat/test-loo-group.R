dense_ll <- function(K, r) {
  U <- chol(K)
  z <- backsolve(U, r, transpose = TRUE)
  -0.5 * sum(z^2) - sum(log(diag(U))) - 0.5 * length(r) * log(2 * pi)
}

test_that("pruning with missing tips equals the marginal likelihood of observed tips", {
  set.seed(41)
  tr <- ape::rtree(20)
  ti <- tree_info(tr)
  k <- k_lambda() + k_white()
  kp <- phylogp:::kernel_parts(k, ti, list(lambda1 = list(lambda = 0.6)))
  K <- kernel_matrix(k, tr, list(lambda1 = list(lambda = 0.6)))
  r <- rnorm(20)
  expect_equal(phylogp:::prune_loglik_na(ti, r, kp$w, kp$rv, kp$t2),
               phylogp:::prune_loglik(ti, r, kp$w, kp$rv, kp$t2))
  miss <- c(2, 7, 8, 15)
  r2 <- r; r2[miss] <- NA
  o <- setdiff(1:20, miss)
  expect_equal(phylogp:::prune_loglik_na(ti, r2, kp$w, kp$rv, kp$t2),
               dense_ll(unname(K[o, o]), r[o]), tolerance = 1e-10)
})

test_that("joint group densities match dense conditioning, with a random root", {
  set.seed(42)
  tr <- ape::rcoal(18)
  ti <- tree_info(tr)
  k <- k_ou() + k_white()
  p <- list(ou1 = list(alpha = 2), white1 = list(sigma = 0.3))
  kp <- phylogp:::kernel_parts(k, ti, p)
  K <- unname(kernel_matrix(k, tr, p))
  r <- rnorm(18)
  idx <- c(3, 4, 5, 11)
  r_rest <- r; r_rest[idx] <- NA
  o <- setdiff(1:18, idx)
  A <- K[idx, o] %*% solve(K[o, o])
  S <- K[idx, idx] - A %*% K[o, idx]
  expect_equal(phylogp:::group_lpd(ti, kp, r, r_rest, kp$t2, idx),
               dense_ll(S, r[idx] - as.numeric(A %*% r[o])), tolerance = 1e-8)
})

test_that("loo_group() refits per clade and works with loo_compare()", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)))
  set.seed(43)
  tr <- ape::rcoal(24)
  y <- kernel_simulate(k_bm() + k_white(), tr)[, 1]
  d <- data.frame(y = y, row.names = names(y))
  clades <- stats::cutree(stats::hclust(stats::as.dist(ape::cophenetic.phylo(tr))), k = 3)
  f1 <- phylo_gp(y ~ 1, d, tr, k_bm() + k_white(), chains = 2, iter_warmup = 200, iter_sampling = 200)
  f2 <- phylo_gp(y ~ 1, d, tr, k_white(), chains = 2, iter_warmup = 200, iter_sampling = 200)
  g1 <- suppressMessages(loo_group(f1, clades))
  g2 <- suppressMessages(loo_group(f2, clades))
  expect_s3_class(g1, "kfold")
  expect_equal(nrow(attr(g1, "groups")), 3)
  expect_equal(unname(g1$estimates["elpd_kfold", "Estimate"]), sum(attr(g1, "groups")$elpd))
  expect_true(all(is.finite(attr(g1, "species_elpd"))))
  cmp <- loo::loo_compare(list(bm = g1, white = g2))
  expect_equal(nrow(cmp), 2)
})
