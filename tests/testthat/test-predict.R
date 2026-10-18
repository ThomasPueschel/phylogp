# Linear-time node posteriors against dense Gaussian conditioning

node_incidence <- function(ti) {
  B <- matrix(0, ti$n_nodes, ti$E)
  parent_edge <- integer(ti$n_nodes)
  parent_edge[ti$des] <- seq_len(ti$E)
  for (n in seq_len(ti$n_nodes)) {
    node <- n
    while (node != ti$root) {
      e <- parent_edge[node]
      B[n, e] <- 1
      node <- ti$anc[e]
    }
  }
  B
}

dense_posterior <- function(C, r, t2) {
  o <- which(!is.na(r))
  Koo <- C[o, o] + diag(t2[o], length(o))
  A <- C[, o] %*% solve(Koo)
  list(mean = as.numeric(A %*% r[o]), var = unname(diag(C) - rowSums(A * C[, o])))
}

test_that("node posteriors match dense conditioning for branch-wise kernels", {
  set.seed(1)
  tr <- ape::rtree(25)                       # not ultrametric
  ti <- tree_info(tr)
  w <- runif(ti$E, 0.1, 1) * ti$len
  for (rv in c(0, 0.7)) {
    t2 <- runif(ti$N, 0.05, 0.3)
    r <- rnorm(ti$N)
    r[sample(ti$N, 6)] <- NA                 # unobserved tips
    B <- node_incidence(ti)
    C <- B %*% (w * t(B)) + rv
    d <- dense_posterior(C, r, c(t2, rep(0, ti$n_nodes - ti$N)))
    p <- phylogp:::tree_posterior(ti, r, w, rv, t2)
    expect_equal(p$mean, d$mean, tolerance = 1e-8)
    expect_equal(p$var, d$var, tolerance = 1e-8)
  }
})

test_that("tips observed without noise are handled exactly", {
  set.seed(2)
  tr <- ape::rcoal(15)
  ti <- tree_info(tr)
  r <- rnorm(ti$N); r[3] <- NA
  B <- node_incidence(ti)
  C <- B %*% (ti$len * t(B))
  d <- dense_posterior(C, r, rep(0, ti$n_nodes))
  p <- phylogp:::tree_posterior(ti, r, ti$len, 0, rep(0, ti$N))
  expect_equal(p$mean, d$mean, tolerance = 1e-8)
  expect_equal(p$var, pmax(d$var, 0), tolerance = 1e-8)
})

test_that("OU transitions reproduce the OU covariance at all nodes", {
  set.seed(3)
  tr <- ape::rtree(20)
  ti <- tree_info(tr)
  h <- ape::node.depth.edgelength(ti$tree)
  Dn <- ape::dist.nodes(ti$tree)
  tn <- (outer(h, h, "+") - Dn) / 2           # height of the MRCA of two nodes
  alpha <- 1.7; s2 <- 2
  a <- exp(-alpha * ti$len)
  w <- s2 / (2 * alpha) * (1 - a^2)
  r <- rnorm(ti$N); r[c(2, 9)] <- NA
  t2 <- rep(0.1, ti$N)
  for (root in c("random", "fixed")) {
    C <- s2 / (2 * alpha) * exp(-alpha * Dn)
    if (root == "fixed") C <- C * (1 - exp(-2 * alpha * tn))
    rv <- if (root == "random") s2 / (2 * alpha) else 0
    d <- dense_posterior(C, r, c(t2, rep(0, ti$n_nodes - ti$N)))
    p <- phylogp:::tree_posterior(ti, r, w, rv, t2, a = a)
    keep <- if (root == "fixed") -ti$root else seq_len(ti$n_nodes)
    expect_equal(p$mean[keep], d$mean[keep], tolerance = 1e-8)
    expect_equal(p$var[keep], d$var[keep], tolerance = 1e-8)
  }
})
