set.seed(1)
tr <- ape::rcoal(30)
ti <- tree_info(tr)
C <- unname(ape::vcv(tr)[ti$labels, ti$labels])
D <- unname(ape::cophenetic.phylo(tr)[ti$labels, ti$labels])
Tm <- ti$Tmax
km <- function(k, pars = NULL) unname(kernel_matrix(k, tr, pars))

test_that("branch-wise kernels equal their closed forms", {
  expect_equal(km(k_bm()), C)
  lam <- 0.3 * C; diag(lam) <- diag(C)
  expect_equal(km(k_lambda(), list(lambda1 = list(lambda = 0.3))), lam)
  expect_equal(km(k_delta(), list(delta1 = list(delta = 2))), C^2 / Tm)
  expect_equal(km(k_eb(), list(eb1 = list(rate = -2))), (exp(-2 * C) - 1) / -2)
  a <- 1.3
  expect_equal(km(k_ou(), list(ou1 = list(alpha = a))), exp(-a * D) / (2 * a))
  expect_equal(km(k_ou(root = "fixed"), list(ou1 = list(alpha = a))),
               exp(-a * D) * (1 - exp(-2 * a * C)) / (2 * a))
})

test_that("OU with a fixed root equals ACDC with rate 2 alpha on an ultrametric tree", {
  a <- 0.8
  ou <- km(k_ou(root = "fixed"), list(ou1 = list(alpha = a)))
  eb <- km(k_eb(upper = 10), list(eb1 = list(sigma = exp(-a * Tm), rate = 2 * a)))
  expect_equal(ou, eb)
})

test_that("lambda equals BM plus white noise on an ultrametric tree", {
  lam <- km(k_lambda(), list(lambda1 = list(lambda = 0.4)))
  bw <- km(k_bm() + k_white(), list(bm1 = list(sigma = sqrt(0.4)),
                                    white1 = list(sigma = sqrt(0.6 * Tm))))
  expect_equal(lam, bw)
})

test_that("user-defined edge kernels reproduce built-in ones", {
  ke <- k_edge("pow(len, kappa)", pars = list(kappa = c(0, 1)))
  expect_equal(km(ke, list(edge1 = list(kappa = 0.4))),
               km(k_kappa(), list(kappa1 = list(kappa = 0.4))))
})

test_that("kernels add and get unique ids", {
  k <- k_bm() + k_white() + k_bm()
  expect_s3_class(k, "phylogp_kernel")
  expect_equal(phylogp:::component_ids(k), c("bm1", "white1", "bm2"))
  expect_equal(km(k), 2 * C + diag(30))
  expect_error(kernel_matrix(k, tr, list(foo1 = list())), "Unknown component")
})

test_that("k_edge rejects reserved names and bad bounds", {
  expect_error(k_edge("len", pars = list(len = c(0, 1))), "Reserved")
  expect_error(k_edge("len", pars = list(a = c(1, 0))), "Bounds")
})
