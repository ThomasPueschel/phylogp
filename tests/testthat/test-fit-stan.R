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
    sd <- kernel_standata(s[[1]], tr, d$y)
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

test_that("loo_refit() agrees with PSIS where PSIS is reliable", {
  skip_if_no_cmdstan()
  set.seed(10)
  tr <- ape::rcoal(30)
  k <- k_bm() + k_white()
  y <- kernel_simulate(k, tr, list(white1 = list(sigma = 0.5)))[, 1]
  d <- data.frame(y = y, row.names = names(y))
  fit <- phylo_gp(y ~ 1, d, tr, k, chains = 2, iter_warmup = 300, iter_sampling = 500)
  lo <- suppressWarnings(loo(fit))
  fake <- lo
  fake$diagnostics$pareto_k[c(3, 17)] <- 1    # pretend two species are influential
  lr <- suppressMessages(loo_refit(fit, loo = fake))
  expect_equal(lr$pointwise[c(3, 17), "elpd_loo"], lo$pointwise[c(3, 17), "elpd_loo"],
               tolerance = 0.05)
  expect_equal(attr(lr, "refitted"), fit$tree$tip.label[c(3, 17)])
  expect_equal(unname(lr$estimates["elpd_loo", "Estimate"]), sum(lr$pointwise[, "elpd_loo"]))
})

# Brute-force check of a left-out species' density at one posterior draw:
# integrate its own latent value numerically against its conditional prior.
loo_by_integration <- function(fit, s, i, own) {
  ti <- tree_info(fit$tree)
  L <- phylogp:::kernel_layout(fit$kernel, ti)
  row <- phylogp:::fit_draws(fit)[s, ]
  kp <- phylogp:::kernel_parts(L$kernel, ti, phylogp:::draw_pars(L, row, ti))
  mu <- as.numeric(fit$standata$X %*% row[grep("^coefs", names(row))])
  r <- fit$standata$y - mu - kp$shift
  r_i <- r[i] + kp$shift[i]                        # residual before i's own shift
  r[i] <- NA
  post <- phylogp:::tree_posterior(ti, r, kp$w, kp$rv, kp$t2)
  log(own(ti, row, kp, r_i, post$mean[i], post$var[i]))
}

test_that("leave-one-out for latent slopes integrates the species' own slope", {
  skip_if_no_cmdstan()
  set.seed(12)
  tr <- ape::rcoal(30)
  x <- kernel_simulate(k_bm(), tr)[, 1]
  y <- 0.5 * x + kernel_simulate(k_bm() + k_white(), tr)[, 1]
  d <- data.frame(y = y, x = x, row.names = names(y))
  fit <- phylo_gp(y ~ x, d, tr, k_bm() + k_slope("x") * k_bm() + k_white(),
                  chains = 1, iter_warmup = 150, iter_sampling = 50)
  ti <- tree_info(fit$tree)
  xs <- fit$standata$slopebm1_xs
  i <- 4; e <- ti$tip_edge[i]
  own <- function(ti, row, kp, r_i, m, v) {
    vp <- row[[paste0("slopebm1_v[", ti$anc[e], "]")]]
    sd_v <- row[["slopebm1_sigma"]] * sqrt(ti$len[e])
    t2 <- kp$t2[i]
    lat_var <- v - kp$w[e] + kp$w[e]              # latent variance of x_i given others
    stats::integrate(function(vi) sapply(vi, function(z)
      stats::dnorm(r_i - xs[i] * z, m, sqrt(lat_var + t2)) * stats::dnorm(z, vp, sd_v)),
      vp - 10 * sd_v, vp + 10 * sd_v, rel.tol = 1e-10)$value
  }
  stan_ll <- fit$fit$draws("log_lik", format = "matrix")[7, i]
  expect_equal(loo_by_integration(fit, 7, i, own), as.numeric(stan_ll), tolerance = 1e-6)
})

test_that("leave-one-out for Student-t residuals integrates the species' own variance", {
  skip_if_no_cmdstan()
  set.seed(13)
  tr <- ape::rcoal(30)
  y <- kernel_simulate(k_bm() + k_white(), tr)[, 1]
  d <- data.frame(y = y, row.names = names(y))
  fit <- phylo_gp(y ~ 1, d, tr, k_bm() + k_white(df = 4),
                  chains = 1, iter_warmup = 150, iter_sampling = 50)
  i <- 9
  own <- function(ti, row, kp, r_i, m, v) {
    s2 <- row[["white1_sigma"]]^2
    dens <- function(lam) stats::dnorm(r_i, m, sqrt(v + s2 * lam)) *
      stats::dgamma(1 / lam, 2, 2) / lam^2         # inverse-gamma(2, 2) density
    stats::integrate(function(l) sapply(l, dens), 0, Inf, rel.tol = 1e-10)$value
  }
  stan_ll <- fit$fit$draws("log_lik", format = "matrix")[5, i]
  expect_equal(loo_by_integration(fit, 5, i, own), as.numeric(stan_ll), tolerance = 1e-5)
})

test_that("dense slope kernels give the R likelihood at the ML estimate", {
  skip_if_no_cmdstan()
  set.seed(14)
  tr <- ape::rcoal(40)
  x <- kernel_simulate(k_bm(), tr)[, 1]
  y <- 1 + 0.5 * x + kernel_simulate(k_bm() + k_white(), tr)[, 1]
  d <- data.frame(y = y, x = x, row.names = names(y))
  for (k in list(k_bm() + k_slope("x", latent = FALSE) * k_bm() + k_white(),
                 k_bm() + k_slope("x") * k_white() + k_white())) {
    fit <- phylo_gp(y ~ x, d, tr, k, method = "optimize")
    est <- fit$fit$mle()
    ids <- phylogp:::component_ids(fit$kernel)
    pars <- stats::setNames(lapply(ids, function(id) list(sigma = est[[paste0(id, "_sigma")]])), ids)
    mu <- est[["coefs[1]"]] + est[["coefs[2]"]] * x
    expect_equal(kernel_loglik(fit$kernel, fit$tree, y, pars, mean = mu[fit$tree$tip.label]),
                 as.numeric(logLik(fit)), tolerance = 1e-6)
  }
})
