dense_loglik <- function(K, y) {
  U <- chol(K)
  z <- backsolve(U, y, transpose = TRUE)
  -0.5 * sum(z^2) - sum(log(diag(U))) - 0.5 * length(y) * log(2 * pi)
}
dense_loo <- function(K, y) {
  Ki <- solve(K)
  g <- as.numeric(Ki %*% y)
  s2 <- 1 / diag(Ki)
  stats::dnorm(y, y - g * s2, sqrt(s2), log = TRUE)
}

kernels <- list(
  lambda_white = k_lambda() + k_white(),
  ou_random = k_ou(),                 # non-zero root variance
  eb = k_eb(),
  bm_dist = k_bm() + k_dist(phi = "exponential") + k_white()  # dense
)

for (nm in names(kernels)) {
  test_that(paste("pruning likelihood and LOO match dense algebra:", nm), {
    set.seed(2)
    tr <- ape::rcoal(25)
    ti <- tree_info(tr)
    k <- kernels[[nm]]
    K <- kernel_matrix(k, tr)
    y <- kernel_simulate(k, tr)[, 1]
    expect_equal(kernel_loglik(k, tr, y), dense_loglik(K, y[ti$labels]))
    kp <- phylogp:::kernel_parts(k, ti)
    if (is.null(kp$K)) {
      expect_equal(phylogp:::prune_loo_lpd(ti, y[ti$labels], kp$w, kp$rv, kp$t2),
                   dense_loo(K, y[ti$labels]), tolerance = 1e-8)
    }
  })
}

test_that("pruning works on non-ultrametric trees with polytomies", {
  set.seed(3)
  tr <- ape::rtree(30)
  tr$edge.length[sample(which(tr$edge[, 2] > 30), 3)] <- 0
  tr <- ape::di2multi(tr)
  ti <- tree_info(tr)
  k <- k_delta() + k_white()
  y <- kernel_simulate(k, tr)[, 1]
  K <- kernel_matrix(k, tr)
  expect_equal(kernel_loglik(k, tr, y), dense_loglik(K, y[ti$labels]))
  kp <- phylogp:::kernel_parts(k, ti)
  expect_equal(phylogp:::prune_loo_lpd(ti, y[ti$labels], kp$w, kp$rv, kp$t2),
               dense_loo(K, y[ti$labels]), tolerance = 1e-8)
})
