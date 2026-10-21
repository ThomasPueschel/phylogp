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

test_that("rows are reduced to species means and deviations", {
  s <- sim_repeated()
  md <- match_data(y ~ x, s$data, s$tree, "species")
  rp <- md$repeated
  expect_false(is.null(rp))
  expect_equal(sum(rp$n), nrow(s$data))
  sp <- md$tree$tip.label
  expect_equal(md$y, as.numeric(tapply(s$data$y, s$data$species, mean)[sp]))
  expect_equal(unname(md$X[, "x"]), as.numeric(tapply(s$data$x, s$data$species, mean)[sp]))
  expect_equal(as.numeric(rowsum(rp$y_dev, rp$species)), rep(0, length(sp)), tolerance = 1e-10)
  expect_null(match_data(y ~ x, s$data[!duplicated(s$data$species), ], s$tree, "species")$repeated)
})

test_that("repeated rows need `species` and cannot have known measurement error", {
  s <- sim_repeated()
  rownames(s$data) <- NULL
  expect_error(match_data(y ~ x, s$data, s$tree, NULL), "species")
  s$data$se <- 0.1
  expect_error(match_data(y ~ x, s$data, s$tree, "species", se = "se"), "cannot be combined")
})

test_that("the Stan program for repeated rows has the within-species part", {
  code <- kernel_stancode(k_bm() + k_white(), repeated = TRUE)
  expect_match(code, "within_logdens(", fixed = TRUE)
  expect_match(code, "real<lower=0> within_sigma;", fixed = TRUE)
  expect_match(code, "t2 += square(within_sigma) ./ n_obs;", fixed = TRUE)
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
