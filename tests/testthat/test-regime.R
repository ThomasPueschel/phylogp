regime_tree <- function(n = 40, seed = 1) {
  set.seed(seed)
  tr <- phytools::pbtree(n = n, scale = 1)
  Q <- matrix(1, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  diag(Q) <- -2
  phytools::sim.history(tr, Q, message = FALSE)
}

# Clade of each edge, as a key shared by trees with different edge orders
edge_keys <- function(t) {
  n <- ape::Ntip(t)
  vapply(t$edge[, 2], function(node)
    paste(sort(t$tip.label[intersect(phytools::getDescendants(t, node), seq_len(n))]),
          collapse = ","), "")
}

test_that("k_regime() is a sum of rate-weighted trees of regime times", {
  skip_if_not_installed("phytools")
  sm <- regime_tree()
  rate <- c(a = 0.2, b = 1, c = 4)
  K <- kernel_matrix(k_regime(reference = "b"), sm,
                     list(regime1 = list(sigma = 1.3, rate_a = 0.2, rate_c = 4)))
  K2 <- Reduce(`+`, lapply(names(rate), function(s) {
    t2 <- sm
    t2$edge.length <- sm$mapped.edge[, s]
    class(t2) <- "phylo"
    1.3^2 * rate[[s]] * ape::vcv.phylo(t2)
  }))
  expect_equal(unname(K), unname(K2[rownames(K), colnames(K)]))
  expect_true(all(kernel_check(k_regime() + k_white(), sm)$valid))
})

test_that("tree_info() and pruning keep the regime map aligned", {
  skip_if_not_installed("phytools")
  sm <- regime_tree()
  ti <- tree_info(sm)
  m <- match(paste(ti$anc, ti$des), paste(sm$edge[, 1], sm$edge[, 2]))
  expect_equal(ti$regime, unname_rows(sm$mapped.edge[m, ]))
  set.seed(2)
  for (i in 1:3) {
    drop <- sample(sm$tip.label, 15)
    p1 <- drop_tips(sm, drop)
    p2 <- phytools::drop.tip.simmap(sm, drop)
    k <- match(edge_keys(p1), edge_keys(p2))
    expect_false(anyNA(k))
    expect_equal(unname(p1$mapped.edge),
                 unname(p2$mapped.edge[k, colnames(p1$mapped.edge)]))
  }
})

test_that("k_regime() reads its regimes from the tree and checks them", {
  skip_if_not_installed("phytools")
  sm <- regime_tree()
  L <- kernel_layout(k_regime(reference = "c"), tree_info(sm))
  expect_equal(L$kernel$components[[1]]$states, c("c", "a", "b"))
  expect_equal(names(L$kernel$components[[1]]$pars), c("rate_a", "rate_b"))
  expect_error(kernel_matrix(k_regime(), ape::rcoal(10)), "mapped regimes")
  expect_error(kernel_matrix(k_regime(states = c("a", "b")), sm), "does not know")
  expect_equal(regime_names(c("open forest", "1", "open-forest")),
               c("open_forest", "1", "open_forest_1"))
})

test_that("maximum likelihood reproduces phytools::brownie.lite()", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "CmdStan not installed")
  skip_if_not_installed("phytools")
  set.seed(3)
  tr <- phytools::pbtree(n = 120, scale = 1)
  Q <- matrix(c(-1, 1, 1, -1), 2, 2, dimnames = list(c("land", "water"), c("land", "water")))
  sm <- phytools::sim.history(tr, Q, message = FALSE)
  y <- phytools::sim.rates(sm, c(land = 1, water = 6))
  f <- phylo_gp(y ~ 1, data.frame(y = y), sm, k_regime(), method = "optimize")
  b <- phytools::brownie.lite(sm, y)
  s2 <- coef(f)[["regime1_sigma"]]^2
  expect_equal(f$fit$lp(), b$logL.multiple, tolerance = 1e-5)
  expect_equal(c(s2, s2 * coef(f)[["regime1_rate_water"]]), unname(b$sig2.multiple),
               tolerance = 1e-3)
})
