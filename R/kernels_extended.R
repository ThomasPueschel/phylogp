# Kernel extensions: relaxed (autocorrelated) rates, rate through time,
# distance-based kernels and user-defined branch-wise kernels.

#' Relaxed (autocorrelated) rates
#'
#' The log rate of evolution itself evolves by Brownian motion along the tree
#' (as in autocorrelated relaxed clocks). The rate of edge e is
#' `sigma^2 * exp((z_start + z_end) / 2)`, with `z` a Brownian motion on the
#' tree with variance `tau^2` per unit time and `z = 0` at the root.
#' `tau = 0` gives constant-rate Brownian motion. The component keeps the
#' branch-wise form, so the likelihood stays linear in the number of tips.
#'
#' The node values of `z` are latent variables, so models with this component
#' can only be fitted by sampling. In leave-one-out densities, the log rate at
#' the left-out tip is integrated out against its conditional prior by
#' Gauss-Hermite quadrature.
#'
#' @param tau_scale Scale of the half-normal prior on `tau`.
#' @return A kernel component.
#' @examples
#' k_relaxed() + k_white()
#' @export
k_relaxed <- function(tau_scale = 1) {
  new_component(
    name = "relaxed", label = "relaxed()", latent = TRUE, tip_latent = TRUE,
    stan_data = "real<lower=0> ID_tau_scale;",
    stan_pars = c("real<lower=0> ID_tau;", "vector[E] ID_eta;"),
    stan_tpars = c(
      "vector[n_nodes] ID_z;",
      "ID_z[root] = 0;",
      "for (k in 1:E) {",
      "  int e = E - k + 1;",
      "  ID_z[des[e]] = ID_z[anc[e]] + ID_tau * sqrt(len[e]) * ID_eta[e];",
      "}"),
    edge = "exp(0.5 * (ID_z[anc] + ID_z[des])) .* len",
    stan_latent_prior = "ID_eta ~ std_normal();",
    stan_prior = "ID_tau ~ normal(0, ID_tau_scale);",
    # tip-edge variance as a function of the tip's own log rate `zi` (edge e)
    stan_tip_loo = "square(ID_sigma) * exp(0.5 * (ID_z[anc[e]] + zi)) * len[e]",
    tip_z = c(parent = "ID_z[anc[e]]", sd = "ID_tau * sqrt(len[e])"),
    standata = function(ti, id, cmp) stats::setNames(list(cmp$tau_scale),
                                                     paste0(id, "_tau_scale")),
    r_eval = function(ti, p, form = "edge") {
      z <- p$z
      if (is.null(z)) z <- simulate_tree_bm(ti, p$tau)
      list(w = p$sigma^2 * exp(0.5 * (z[ti$anc] + z[ti$des])) * ti$len, rv = 0)
    },
    defaults = function(ti) list(sigma = 1, tau = 1),
    validity = function(ti, form) list(guaranteed = TRUE,
      reason = "branch-wise form; edge variances are positive for every rate path"),
    extra = list(tau_scale = tau_scale)
  )
}

#' Rate of evolution as a smooth function of time
#'
#' The rate `sigma^2(t)` is piecewise constant over `n_bins` equal time bins
#' from the root to the tips, with a Gaussian random-walk prior on its
#' logarithm. Each edge accumulates the integral of the rate over its time
#' span, so the component keeps the branch-wise form. Early burst is a special
#' case (a linear trend in log rate).
#'
#' The bin values are latent variables, so models with this component can
#' only be fitted by sampling.
#'
#' @param n_bins Number of time bins.
#' @param tau_scale Scale of the half-normal prior on the random-walk standard
#'   deviation `tau`.
#' @return A kernel component.
#' @examples
#' k_rate_time(n_bins = 20) + k_white()
#' @export
k_rate_time <- function(n_bins = 20, tau_scale = 0.5) {
  n_bins <- as.integer(n_bins)
  if (n_bins < 2) stop("`n_bins` must be at least 2.", call. = FALSE)
  new_component(
    name = "ratetime", label = sprintf("rate_time(n_bins = %d)", n_bins),
    latent = TRUE,
    stan_data = c("int<lower=2> ID_M;", "matrix<lower=0>[E, ID_M] ID_O;",
                  "real<lower=0> ID_tau_scale;"),
    stan_pars = c("real<lower=0> ID_tau;", "vector[ID_M - 1] ID_eps;"),
    stan_tpars = c(
      "vector[ID_M] ID_lograte;",
      "ID_lograte[1] = 0;",
      "for (m in 2:ID_M) ID_lograte[m] = ID_lograte[m - 1] + ID_tau * ID_eps[m - 1];"),
    edge = "ID_O * exp(ID_lograte)",
    stan_latent_prior = "ID_eps ~ std_normal();",
    stan_prior = "ID_tau ~ normal(0, ID_tau_scale);",
    standata = function(ti, id, cmp) stats::setNames(
      list(cmp$n_bins, bin_overlap(ti, cmp$n_bins), cmp$tau_scale),
      paste0(id, c("_M", "_O", "_tau_scale"))),
    r_eval = function(ti, p, form = "edge") {
      f <- p$lograte
      if (is.null(f)) f <- c(0, cumsum(p$tau * stats::rnorm(n_bins - 1)))
      O <- bin_overlap(ti, length(f))
      list(w = p$sigma^2 * as.numeric(O %*% exp(f)), rv = 0)
    },
    defaults = function(ti) list(sigma = 1, tau = 0.5),
    validity = function(ti, form) list(guaranteed = TRUE,
      reason = "branch-wise form; edge variances are integrals of a positive rate"),
    extra = list(n_bins = n_bins, tau_scale = tau_scale)
  )
}

# Time each edge spends in each of M equal bins between the root and Tmax
bin_overlap <- function(ti, M) {
  brk <- seq(0, ti$Tmax, length.out = M + 1)
  O <- vapply(seq_len(M), function(m)
    pmax(0, pmin(ti$he, brk[m + 1]) - pmax(ti$hs, brk[m])), numeric(ti$E))
  matrix(O, ti$E, M)
}

#' Distance-based kernels
#'
#' A kernel `sigma^2 * phi(D / ell)` on a distance matrix: patristic distances
#' on the tree (the default) or any other distance between species, such as
#' great-circle distances between range centroids for a phylogeny + space
#' model (`k_bm() + k_dist(geo) + k_white()`). The component is dense, so
#' models containing it use dense linear algebra (O(N^3)).
#'
#' Validity depends on the distance and on `phi`. On an ultrametric tree every
#' non-increasing `phi` gives a valid kernel on patristic distances. On trees
#' with fossil tips or without a clock only completely monotone functions such
#' as the exponential are guaranteed valid; squared-exponential and Matern
#' kernels can fail. Use [kernel_check()] before fitting.
#'
#' @param D `NULL` for patristic distances on the tree, or a symmetric
#'   distance matrix with row and column names matching the tip labels.
#' @param phi `"exponential"`, `"matern32"`, `"matern52"` or `"sqexp"`.
#' @param lower,upper Bounds for the length-scale `ell`, as multiples of the
#'   largest distance.
#' @return A kernel component.
#' @examples
#' k_dist(phi = "exponential")
#' @export
k_dist <- function(D = NULL, phi = c("exponential", "matern32", "matern52", "sqexp"),
                   lower = 1e-3, upper = 10) {
  phi <- match.arg(phi)
  if (!is.null(D)) {
    D <- as.matrix(D)
    if (is.null(rownames(D)) || !isSymmetric(unname(D)))
      stop("`D` must be a symmetric matrix with tip labels as row names.", call. = FALSE)
  }
  phi_stan <- switch(phi,
    exponential = "exp(-RR)",
    matern32 = "(1 + sqrt(3) * RR) .* exp(-sqrt(3) * RR)",
    matern52 = "(1 + sqrt(5) * RR + 5 * square(RR) / 3) .* exp(-sqrt(5) * RR)",
    sqexp = "exp(-0.5 * square(RR))")
  maxD <- function(ti, cmp) max(cmp_distance(ti, cmp))
  new_component(
    name = "dist", label = sprintf("dist(%s%s)", if (is.null(D)) "" else "D, ", phi),
    pars = list(ell = bounded(function(ti, cmp) lower * maxD(ti, cmp),
                              function(ti, cmp) upper * maxD(ti, cmp))),
    form = function(ti) "dense",
    dense = phi_stan,
    dense_alias = "matrix[N, N] RR = ID_D / ell;",
    stan_data = c("matrix[N, N] ID_D;", "real<lower=0> ID_ell_loc;"),
    stan_prior = "ID_ell ~ lognormal(log(ID_ell_loc), 1);",
    standata = function(ti, id, cmp) {
      Dm <- cmp_distance(ti, cmp)
      stats::setNames(list(Dm, stats::median(Dm[upper.tri(Dm)])),
                      paste0(id, c("_D", "_ell_loc")))
    },
    r_eval = function(ti, p, form = "dense") {
      RR <- cmp_distance(ti, p$cmp) / p$ell
      K <- switch(phi,
        exponential = exp(-RR),
        matern32 = (1 + sqrt(3) * RR) * exp(-sqrt(3) * RR),
        matern52 = (1 + sqrt(5) * RR + 5 * RR^2 / 3) * exp(-sqrt(5) * RR),
        sqexp = exp(-0.5 * RR^2))
      list(K = p$sigma^2 * K)
    },
    defaults = function(ti) list(sigma = 1, ell = NA),
    validity = function(ti, form) {
      user_D <- !is.null(D)
      if (!user_D && ti$ultrametric) {
        list(guaranteed = TRUE,
             reason = "non-increasing function of patristic distance on an ultrametric tree")
      } else if (!user_D && phi == "exponential") {
        list(guaranteed = TRUE,
             reason = "completely monotone function of a tree metric (valid on every tree)")
      } else if (!user_D) {
        list(guaranteed = FALSE,
             reason = sprintf("%s kernel on patristic distances of a non-ultrametric tree is not guaranteed valid", phi))
      } else {
        list(guaranteed = FALSE,
             reason = sprintf("%s kernel on a user-supplied distance: validity depends on the distance (checked numerically)", phi))
      }
    },
    base_var = function(ti) 1,
    extra = list(D = D, phi = phi)
  )
}

# Distance matrix for a k_dist component, in tip order
cmp_distance <- function(ti, cmp) {
  if (is.null(cmp$D)) return(tree_D(ti))
  miss <- setdiff(ti$labels, rownames(cmp$D))
  if (length(miss))
    stop("Distance matrix lacks ", length(miss), " tips, e.g. ", miss[1], call. = FALSE)
  unname(cmp$D[ti$labels, ti$labels])
}

#' User-defined branch-wise kernels
#'
#' Defines a new tree kernel by the variance each edge accumulates, written
#' as a Stan expression. Any kernel of this form is a valid covariance on
#' every tree as long as the edge variances are non-negative, and it keeps
#' the O(N) pruning likelihood and leave-one-out densities.
#'
#' The expression must return a vector with one element per edge and may use
#' the edge vectors `len` (length), `hs` and `he` (heights of the start and
#' end of the edge above the root), `is_tip` (1 for tip edges), the tree
#' height `Tmax`, and the parameters named in `pars`. Use Stan's elementwise
#' operators `.*` and `./` between vectors. The same expression is evaluated
#' in R (after translating `.*`, `./`, `pow()` and `square()`) for
#' [kernel_matrix()], [kernel_check()] and [kernel_plot()].
#'
#' @param w Stan expression for the edge variances (before scaling by
#'   `sigma^2`).
#' @param pars Named list of parameter bounds, each `c(lower, upper)`.
#'   Parameters get uniform priors on these bounds.
#' @param defaults Named list of parameter values used by [kernel_plot()] and
#'   [kernel_check()]; defaults to the midpoints of the bounds.
#' @param name Short name used in printed output and Stan variable names.
#' @return A kernel component.
#' @examples
#' # Pagel's kappa with a rate trend through time
#' k <- k_edge("pow(len, kappa) .* exp(beta * he)",
#'             pars = list(kappa = c(0, 1), beta = c(-5, 5)))
#' k
#' @export
k_edge <- function(w, pars = list(), defaults = NULL, name = "edge") {
  if (!is.character(w) || length(w) != 1) stop("`w` must be one string.", call. = FALSE)
  if (!grepl("^[a-z][a-z0-9]*$", name)) stop("`name` must be lower-case letters and digits.", call. = FALSE)
  reserved <- c("len", "hs", "he", "is_tip", "Tmax", "N", "P", "E", "X", "y", "se",
                "n_nodes", "root", "anc", "des", "first_child", "tip_edge", "B",
                "tC", "Dtree", "use_prior", "prior_beta_scale", "coefs", "w", "rv",
                "t2", "K", "r", "mo", "log_lik", "sigma", "e", "i", "k", "g", "m",
                "zi", "zp", "sdz", "lp", "n_gh", "gh_x", "gh_w")
  bad <- intersect(names(pars), reserved)
  if (length(bad)) stop("Reserved parameter name: ", bad[1], call. = FALSE)
  for (p in names(pars)) {
    b <- pars[[p]]
    if (length(b) != 2 || !all(is.finite(b)) || b[1] >= b[2])
      stop("Bounds for `", p, "` must be c(lower, upper) with finite values.", call. = FALSE)
  }
  r_expr <- stan_to_r(w)
  defaults <- defaults %||% lapply(pars, mean)
  new_component(
    name = name, label = sprintf("edge(\"%s\")", w),
    pars = lapply(pars, function(b) bounded(b[1], b[2])),
    edge = w,
    r_eval = function(ti, p, form = "edge") {
      env <- c(list(len = ti$len, hs = ti$hs, he = ti$he, is_tip = ti$is_tip,
                    Tmax = ti$Tmax), p[names(pars)])
      ww <- eval(r_expr, envir = env, enclos = stan_r_env)
      if (length(ww) == 1) ww <- rep(ww, ti$E)
      list(w = p$sigma^2 * ww, rv = 0)
    },
    defaults = function(ti) c(list(sigma = 1), defaults),
    validity = function(ti, form) list(guaranteed = NA,
      reason = "branch-wise form: valid wherever all edge variances are non-negative (checked numerically)")
  )
}

# Translate a simple Stan expression into R
stan_to_r <- function(x) {
  x <- gsub(".*", "*", x, fixed = TRUE)
  x <- gsub("./", "/", x, fixed = TRUE)
  tryCatch(str2lang(x), error = function(e)
    stop("Could not translate the edge expression to R: ", conditionMessage(e), call. = FALSE))
}

stan_r_env <- local({
  env <- new.env(parent = baseenv())
  env$pow <- function(x, y) x^y
  env$square <- function(x) x^2
  env$cube <- function(x) x^3
  env$inv <- function(x) 1 / x
  env$inv_logit <- stats::plogis
  env$fmin <- pmin
  env$fmax <- pmax
  env$step <- function(x) as.numeric(x > 0)
  env$log1m <- function(x) log1p(-x)
  env
})

# One Brownian path on the tree (z = 0 at the root), variance tau^2 per unit time
simulate_tree_bm <- function(ti, tau) {
  z <- numeric(ti$n_nodes)
  for (e in rev(seq_len(ti$E)))
    z[ti$des[e]] <- z[ti$anc[e]] + tau * sqrt(ti$len[e]) * stats::rnorm(1)
  z
}

`%||%` <- function(x, y) if (is.null(x)) y else x
