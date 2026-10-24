# Several rows per species

sim_repeated <- function(seed = 1, n_sp = 30) {
  set.seed(seed)
  tr <- ape::rcoal(n_sp)
  f <- kernel_simulate(k_bm(), tr)[, 1] + rnorm(n_sp, 0, 0.3)
  n <- sample(1:4, n_sp, replace = TRUE)
  sp <- rep(names(f), n)
  x <- rnorm(length(sp))
  d <- data.frame(species = sp, x = x, y = 1 + 0.5 * x + f[sp] + rnorm(length(sp), 0, 0.5))
  list(tree = tr, data = d[sample(nrow(d)), ])
}

test_that("rows are reduced to species means", {
  s <- sim_repeated()
  md <- match_data(y ~ x, s$data, s$tree, "species")
  rp <- md$repeated
  expect_false(is.null(rp))
  expect_equal(sum(rp$n), nrow(s$data))
  sp <- md$tree$tip.label
  expect_equal(md$y, as.numeric(tapply(s$data$y, s$data$species, mean)[sp]))
  expect_equal(unname(md$X[, "x"]), as.numeric(tapply(s$data$x, s$data$species, mean)[sp]))
  expect_equal(rp$y, rp$rows$y)
  expect_equal(sp[rp$species], as.character(rp$rows$species))
  expect_null(match_data(y ~ x, s$data[!duplicated(s$data$species), ], s$tree, "species")$repeated)
})

test_that("repeated rows need `species` and keep the measurement error of each row", {
  s <- sim_repeated()
  rownames(s$data) <- NULL
  expect_error(match_data(y ~ x, s$data, s$tree, NULL), "species")
  s$data$se <- runif(nrow(s$data), 0.05, 0.2)
  rp <- match_data(y ~ x, s$data, s$tree, "species", "se", se = "se")$repeated
  expect_equal(rp$se, rp$rows$se)
})

test_that("within_parts() factorises the likelihood of the rows exactly", {
  set.seed(4)
  sp <- c(1, 1, 2, 3, 3, 3)
  r <- rnorm(6); se <- runif(6, 0, 0.5); s <- 0.7; f <- c(0.3, -0.2, 0.5)
  wp <- within_parts(r, se, sp, 3, s)
  direct <- as.numeric(rowsum(dnorm(r, f[sp], sqrt(s^2 + se^2), log = TRUE), sp))
  expect_equal(dnorm(wp$rbar, f, sqrt(wp$v), log = TRUE) + wp$lp, direct)
})

test_that("the Stan program for repeated rows has the within-species part", {
  code <- kernel_stancode(k_bm() + k_white(), repeated = TRUE)
  expect_match(code, "within_parts(", fixed = TRUE)
  expect_match(code, "real<lower=0> within_sigma;", fixed = TRUE)
  expect_match(code, "t2 += wp[, 2];", fixed = TRUE)
  expect_false(grepl("within_sigma", kernel_stancode(k_bm() + k_white())))
})

test_that("the likelihood with repeated rows equals the dense likelihood of all rows", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "CmdStan not installed")
  s <- sim_repeated(2, 25)
  d <- s$data
  fit <- suppressMessages(phylo_gp(y ~ x, d, s$tree, k_bm() + k_white(),
                                   species = "species", method = "optimize"))
  m <- fit$fit$mle()
  sp <- fit$tree$tip.label
  C <- ape::vcv.phylo(fit$tree)[sp, sp]
  Z <- outer(d$species, sp, "==") * 1
  S <- m[["bm1_sigma"]]^2 * Z %*% C %*% t(Z) + m[["white1_sigma"]]^2 * tcrossprod(Z) +
    diag(m[["within_sigma"]]^2, nrow(d))
  r <- d$y - cbind(1, d$x) %*% m[c("coefs[1]", "coefs[2]")]
  U <- chol(S)
  z <- backsolve(U, r, transpose = TRUE)
  dense <- -0.5 * sum(z^2) - sum(log(diag(U))) - 0.5 * nrow(d) * log(2 * pi)
  expect_equal(as.numeric(logLik(fit)), dense, tolerance = 1e-6)
  expect_true("within_sigma" %in% summary(fit)$parameter)
})

test_that("with known errors per row, the likelihood equals the dense likelihood", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "CmdStan not installed")
  s <- sim_repeated(3, 25)
  d <- s$data
  d$se <- runif(nrow(d), 0.05, 0.6)
  fit <- suppressMessages(phylo_gp(y ~ x, d, s$tree, k_bm() + k_white(), species = "species",
                                   se = "se", method = "optimize"))
  m <- fit$fit$mle()
  sp <- fit$tree$tip.label
  C <- ape::vcv.phylo(fit$tree)[sp, sp]
  Z <- outer(d$species, sp, "==") * 1
  S <- m[["bm1_sigma"]]^2 * Z %*% C %*% t(Z) + m[["white1_sigma"]]^2 * tcrossprod(Z) +
    diag(m[["within_sigma"]]^2 + d$se^2)
  r <- d$y - cbind(1, d$x) %*% m[c("coefs[1]", "coefs[2]")]
  U <- chol(S)
  z <- backsolve(U, r, transpose = TRUE)
  dense <- -0.5 * sum(z^2) - sum(log(diag(U))) - 0.5 * nrow(d) * log(2 * pi)
  expect_equal(as.numeric(logLik(fit)), dense, tolerance = 1e-6)
})

test_that("within_parts() integrates a species slope exactly", {
  set.seed(5)
  sp <- c(1, 1, 1, 2, 2, 3)
  r <- rnorm(6); se <- runif(6, 0, 0.5); s <- 0.6; xs <- rnorm(6)
  s2u <- c(0.4, 0.9, 0.2); g <- c(0.1, -0.3, 0.2)
  wp <- within_parts(r, se, sp, 3, s, xs, s2u)
  direct <- vapply(1:3, function(i) {
    k <- sp == i
    S <- diag(s^2 + se[k]^2, sum(k)) + s2u[i] * tcrossprod(xs[k])
    U <- chol(S); z <- backsolve(U, r[k] - g[i], transpose = TRUE)
    -0.5 * sum(z^2) - sum(log(diag(U))) - 0.5 * sum(k) * log(2 * pi)
  }, 0)
  expect_equal(dnorm(wp$rbar, g, sqrt(wp$v), log = TRUE) + wp$lp, direct)
})

test_that("varying slopes with repeated rows: leave-one-out matches an R calculation", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
          "CmdStan not installed")
  s <- sim_repeated(6, 20)
  d <- s$data
  d$se <- runif(nrow(d), 0, 0.3)
  fit <- suppressWarnings(suppressMessages(phylo_gp(
    y ~ x, d, s$tree, k_bm() + k_slope("x") * k_bm() + k_white(), species = "species",
    se = "se", chains = 1, iter_warmup = 200, iter_sampling = 20)))
  expect_match(fit$code, "slopebm1_xs_obs", fixed = TRUE)
  ll <- fit$fit$draws("log_lik", format = "matrix")
  D <- fit_draws(fit)
  ti <- tree_info(fit$tree)
  L <- fit$layout
  sd <- fit$standata
  for (k in c(1, 10)) {
    row <- D[k, ]
    beta <- unname(row[c("coefs[1]", "coefs[2]")])
    kp <- kernel_parts(L$kernel, ti, draw_pars(L, row, ti))
    sr <- species_resid(fit, row, beta, kp)
    v <- kp$pars$slopebm1$v
    wb <- L$kernel$components[[2]]$base$r_eval(ti, kp$pars$slopebm1, "edge")$w
    rr <- sd$y_obs - as.numeric(sd$X_obs %*% beta) -
      sd$slopebm1_xs_obs * v[ti$anc[ti$tip_edge]][sd$obs_sp]
    wl <- within_parts(rr, sd$se_obs, sd$obs_sp, sd$N, unname(row["within_sigma"]),
                       sd$slopebm1_xs_obs, wb[ti$tip_edge])
    t2 <- kp$t2 + sr$t2
    ref <- vapply(seq_len(ti$N), function(i) {
      r <- sr$r; r[i] <- NA
      cm <- condition_tips(ti, kp, r, t2, i)
      dnorm(wl$rbar[i], cm$mean, sqrt(cm$var + kp$t2[i] + wl$v[i]), log = TRUE) + wl$lp[i]
    }, 0)
    expect_equal(as.numeric(ll[k, ]), ref, tolerance = 1e-6)
  }
})
