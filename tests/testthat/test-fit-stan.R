# These tests compile Stan programs; they run only when CmdStan is available
# and not on CRAN.

skip_if_no_cmdstan <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "CmdStan not installed")
}

test_that("maximum likelihood reproduces phylolm for the classical models", {
  skip_if_no_cmdstan()
  skip_if_not_installed("phylolm")
  set.seed(3)
  tr <- ape::rcoal(60)
  x <- kernel_simulate(k_bm(), tr)[, 1]
  specs <- list(
    list(k_bm(), "BM", NULL),
    list(k_lambda(), "lambda", list(lambda1 = list(lambda = 0.6))),
    list(k_kappa(), "kappa", list(kappa1 = list(kappa = 0.4))),
    list(k_delta(), "delta", list(delta1 = list(delta = 2))),
    list(k_eb(), "EB", list(eb1 = list(rate = -1))),
    list(k_ou(root = "fixed"), "OUfixedRoot", list(ou1 = list(alpha = 1))),
    list(k_ou(), "OUrandomRoot", list(ou1 = list(alpha = 1))))
  for (s in specs) {
    e <- kernel_simulate(s[[1]], tr, s[[3]])[, 1]
    d <- data.frame(y = 1 + 0.5 * x + e, x = x, row.names = names(e))
    sd <- make_standata(s[[1]], tr, d$y)
    args <- list(y ~ x, data = d, phy = tr, model = s[[2]])
    lb <- grep("_lb$", names(sd))
    if (length(lb)) {
      args$lower.bound <- sd[[lb]]
      args$upper.bound <- sd[[grep("_ub$", names(sd))]]
    }
    pl <- suppressWarnings(do.call(phylolm::phylolm, args))
    fit <- suppressMessages(phylo_gp(y ~ x, d, tr, s[[1]], method = "optimize"))
    expect_equal(as.numeric(logLik(fit)), pl$logLik, tolerance = 1e-5,
                 label = paste("logLik", s[[2]]))
    expect_equal(unname(coef(fit)[2]), unname(coef(pl)[2]), tolerance = 1e-3,
                 label = paste("slope", s[[2]]))
  }
})

test_that("leave-one-out densities from Stan match the R recursion", {
  skip_if_no_cmdstan()
  set.seed(5)
  tr <- ape::rcoal(40)
  ti <- tree_info(tr)
  k <- k_ou() + k_white()
  y <- kernel_simulate(k, tr)[, 1]
  d <- data.frame(y = y, row.names = names(y))
  fit <- phylo_gp(y ~ 1, d, tr, k, chains = 2, iter_warmup = 200, iter_sampling = 100)
  dr <- as.data.frame(posterior::as_draws_df(fit$fit$draws()))[10, ]
  pars <- list(ou1 = list(sigma = dr$ou1_sigma, alpha = dr$ou1_alpha),
               white1 = list(sigma = dr$white1_sigma))
  kp <- phylogp:::kernel_parts(k, ti, pars)
  r <- y[ti$labels] - dr[["coefs[1]"]]
  expect_equal(unname(phylogp:::prune_loo_lpd(ti, r, kp$w, kp$rv, kp$t2)),
               unname(unlist(dr[paste0("log_lik[", 1:40, "]")])), tolerance = 1e-6)
  expect_s3_class(loo(fit), "loo")
})

test_that("ancestral states reproduce phytools::fastAnc under BM", {
  skip_if_no_cmdstan()
  skip_if_not_installed("phytools")
  set.seed(8)
  tr <- ape::rcoal(30)
  y <- kernel_simulate(k_bm(), tr)[, 1] + 2
  d <- data.frame(y = y, row.names = names(y))
  fit <- phylo_gp(y ~ 1, d, tr, k_bm(), method = "optimize")
  anc <- ancestral(fit)
  fa <- phytools::fastAnc(fit$tree, y[fit$tree$tip.label])
  expect_equal(anc$mean, unname(fa[as.character(anc$node)]), tolerance = 1e-5)
})

test_that("predict() for unobserved species equals dense kriging", {
  skip_if_no_cmdstan()
  set.seed(9)
  tr <- ape::rcoal(40)
  k <- k_lambda()
  x <- kernel_simulate(k_bm(), tr)[, 1]
  y <- 1 + 0.5 * x + kernel_simulate(k, tr, list(lambda1 = list(lambda = 0.7)))[, 1]
  d <- data.frame(y = y, x = x, row.names = names(y))
  held <- names(y)[1:5]
  fit <- suppressMessages(phylo_gp(y ~ x, d[setdiff(names(y), held), ], tr, k,
                                   method = "optimize"))
  pr <- predict(fit, newdata = d[held, "x", drop = FALSE])
  est <- fit$fit$mle()
  K <- kernel_matrix(k, tr, list(lambda1 = list(sigma = est[["lambda1_sigma"]],
                                                lambda = est[["lambda1_lambda"]])))
  o <- setdiff(names(y), held)
  mu <- est[["coefs[1]"]] + est[["coefs[2]"]] * x
  A <- K[held, o] %*% solve(K[o, o])
  expect_equal(pr$mean, as.numeric(mu[held] + A %*% (y[o] - mu[o])), tolerance = 1e-6)
  expect_equal(pr$sd, unname(sqrt(diag(K[held, held] - A %*% K[o, held]))), tolerance = 1e-6)
  expect_false(any(pr$observed))
})
