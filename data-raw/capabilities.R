# Which features work with which kernel components: fits every kernel to a
# small simulated data set and tries each feature, recording "ok" or the error.
# Run from the package root with the package installed:
#   Rscript data-raw/capabilities.R
# Writes data-raw/capabilities.csv (used by the capabilities table in the docs).

suppressPackageStartupMessages({ library(ape); library(phylogp) })
set.seed(1)
tr <- rcoal(24)
sp <- tr$tip.label
x <- setNames(rnorm(24), sp)
y <- 1 + kernel_simulate(k_bm(), tr)[, 1] + 0.5 * x + rnorm(24, 0, 0.3)
d <- data.frame(species = sp, y = y[sp], x = x, se = 0.1, row.names = sp)
fit_sp <- sp[-1]                       # sp[1] is in the tree but not fitted
d_fit <- d[fit_sp, ]
# several rows per species
n_rep <- rep(1:3, length.out = length(fit_sp))
d_rep <- d_fit[rep(seq_along(fit_sp), n_rep), ]
d_rep$y <- d_rep$y + rnorm(nrow(d_rep), 0, 0.2)
rownames(d_rep) <- NULL
trees <- c(tr, rcoal(24, tip.label = sample(sp)))
groups <- setNames(rep(1:4, length.out = length(fit_sp)), fit_sp)

kernels <- list(
  "k_bm() + k_white()" = k_bm() + k_white(),
  "k_lambda()" = k_lambda(),
  "k_kappa() + k_white()" = k_kappa() + k_white(),
  "k_delta() + k_white()" = k_delta() + k_white(),
  "k_eb() + k_white()" = k_eb() + k_white(),
  "k_ou() + k_white()" = k_ou() + k_white(),
  "k_bm() + k_white(df = 4)" = k_bm() + k_white(df = 4),
  "k_relaxed() + k_white()" = k_relaxed() + k_white(),
  "k_rate_time(5) + k_white()" = k_rate_time(5) + k_white(),
  "k_bm() + k_dist() + k_white()" = k_bm() + k_dist() + k_white(),
  "k_edge() + k_white()" = k_edge("pow(len, kappa)", pars = list(kappa = c(0, 1))) + k_white(),
  "k_bm() + k_slope(x) * k_bm() + k_white()" = k_bm() + k_slope("x") * k_bm() + k_white(),
  "k_bm() + k_slope(x, latent = FALSE) * k_bm() + k_white()" =
    k_bm() + k_slope("x", latent = FALSE) * k_bm() + k_white(),
  "k_bm() + k_slope(x) * k_white()" = k_bm() + k_slope("x") * k_white())

hmc <- list(chains = 1, iter_warmup = 150, iter_sampling = 150)
try1 <- function(expr) {
  r <- tryCatch({ suppressWarnings(suppressMessages(force(expr))); "ok" },
                error = function(e) conditionMessage(e))
  gsub("\\s+", " ", r)
}
fit1 <- function(data, ...) suppressWarnings(suppressMessages(
  do.call(phylo_gp, c(list(y ~ 1, data, tr, k, species = "species"), list(...), hmc))))

rows <- list()
for (nm in names(kernels)) {
  k <- kernels[[nm]]
  message(nm)
  res <- list(kernel = nm)
  res$ml <- try1(suppressWarnings(suppressMessages(
    phylo_gp(y ~ 1, d_fit, tr, k, species = "species", method = "optimize"))))
  f <- tryCatch(fit1(d_fit), error = function(e) e)
  res$hmc <- if (inherits(f, "error")) conditionMessage(f) else "ok"
  ok <- !inherits(f, "error")
  res$loo <- if (ok) try1(loo(f)) else NA
  # refit the species with the largest Pareto k
  res$loo_refit <- if (ok) try1({
    lo <- loo(f)
    loo_refit(f, lo, k_threshold = max(lo$diagnostics$pareto_k) - 1e-9)
  }) else NA
  res$loo_group <- if (ok) try1(loo_group(f, groups, ndraws = 50)) else NA
  res$predict_fitted <- if (ok) try1(predict(f, ndraws = 50)) else NA
  res$predict_new <- if (ok) try1(predict(f, newdata = d[sp[1], ], ndraws = 50)) else NA
  res$ancestral <- if (ok) try1(ancestral(f, ndraws = 50)) else NA
  res$se <- try1(fit1(d_fit, se = "se"))
  fr <- tryCatch(fit1(d_rep), error = function(e) e)
  res$repeated <- if (inherits(fr, "error")) conditionMessage(fr) else "ok"
  okr <- !inherits(fr, "error")
  res$repeated_loo <- if (okr) try1(loo(fr)) else NA
  res$repeated_predict_new <- if (okr) try1(predict(fr, newdata = d[sp[1], ], ndraws = 50)) else NA
  res$repeated_se <- try1(fit1(transform(d_rep, se = 0.1), se = "se"))
  res$trees <- try1({
    ft <- do.call(phylo_gp_trees, c(list(y ~ 1, d_fit, trees, k, species = "species"), hmc))
    loo(ft)
  })
  rows[[nm]] <- as.data.frame(res, check.names = FALSE)
}
out <- do.call(rbind, rows)
rownames(out) <- NULL
write.csv(out, "data-raw/capabilities.csv", row.names = FALSE)
print(out[, -1] == "ok")
