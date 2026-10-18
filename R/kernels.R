# Kernel components and the kernel grammar.
#
# A kernel is a sum of components. Each component knows
#   * how it enters the tree covariance: through edge variances w, the root
#     variance rv and extra tip variance t2 (the branch-wise form, which keeps
#     the O(N) pruning likelihood), or as a dense N x N matrix;
#   * the Stan code that declares its parameters and adds its contribution;
#   * an R evaluator, used by kernel_matrix(), kernel_check() and kernel_plot().
# Every component has its own scale parameter `sigma` (its contribution is
# sigma^2 times a base kernel).

new_component <- function(name, label, pars = list(), form = function(ti) "edge",
                          edge = NULL, root = NULL, dense = NULL, white = FALSE,
                          latent = FALSE, tip_latent = FALSE, needs = character(),
                          stan_data = NULL, stan_pars = NULL, stan_tpars = NULL,
                          stan_latent_prior = NULL, stan_prior = NULL,
                          stan_tip_loo = NULL, standata = NULL,
                          tip_z = NULL, dense_alias = NULL, white_expr = NULL,
                          stan_shift = NULL, tip_loo = NULL, hyper = character(),
                          draw_pars = NULL,
                          r_eval = NULL, defaults = NULL, validity = NULL,
                          base_var = function(ti) ti$Tmax, extra = list()) {
  # tip_loo: how leave-one-out treats the left-out species' own latent values:
  #   list(kind = "closed", r = <added to r_i>, v = <added to variance>)
  #   list(kind = "gh" | "laguerre", cur = <current tip variance>,
  #        node = <tip variance at quadrature node g>, setup = <code>)
  structure(c(list(
    name = name, label = label, pars = pars, form = form,
    edge = edge, root = root, dense = dense, white = white,
    latent = latent, tip_latent = tip_latent, needs = needs,
    stan_data = stan_data, stan_pars = stan_pars, stan_tpars = stan_tpars,
    stan_latent_prior = stan_latent_prior, stan_prior = stan_prior,
    stan_tip_loo = stan_tip_loo, standata = standata,
    tip_z = tip_z, dense_alias = dense_alias, white_expr = white_expr,
    stan_shift = stan_shift, tip_loo = tip_loo, hyper = hyper, draw_pars = draw_pars,
    r_eval = r_eval, defaults = defaults, validity = validity,
    base_var = base_var
  ), extra), class = "phylogp_component")
}

as_kernel <- function(x) {
  if (inherits(x, "phylogp_kernel")) return(x)
  if (inherits(x, "phylogp_component")) {
    return(structure(list(components = list(x)), class = "phylogp_kernel"))
  }
  stop("Not a phylogp kernel.", call. = FALSE)
}

#' Combine kernels
#'
#' Adding kernels describes independent processes acting together; the
#' covariance of the sum is the sum of the covariances.
#'
#' @param e1,e2 Kernels or kernel components.
#' @return A kernel of class `phylogp_kernel`.
#' @examples
#' k <- k_bm() + k_white()    # Brownian motion plus residual variance
#' k
#' @name kernel-sum
#' @export
"+.phylogp_kernel" <- function(e1, e2) {
  if (missing(e2)) return(e1)
  structure(list(components = c(as_kernel(e1)$components, as_kernel(e2)$components)),
            class = "phylogp_kernel")
}

#' @rdname kernel-sum
#' @export
"+.phylogp_component" <- `+.phylogp_kernel`

# Unique ids (bm1, white1, ...) used as prefixes for Stan names
component_ids <- function(kernel) {
  nm <- vapply(kernel$components, `[[`, "", "name")
  ave_idx <- stats::ave(seq_along(nm), nm, FUN = seq_along)
  paste0(nm, ave_idx)
}

#' @export
print.phylogp_kernel <- function(x, ...) {
  ids <- component_ids(x)
  cat("phylogp kernel:",
      paste(vapply(x$components, `[[`, "", "label"), collapse = " + "), "\n")
  for (j in seq_along(x$components)) {
    cmp <- x$components[[j]]
    pn <- c("sigma", names(cmp$pars), cmp$hyper, if (cmp$latent) "(latent values)")
    cat(sprintf("  %-10s %-28s parameters: %s\n", ids[j], cmp$label,
                paste(pn, collapse = ", ")))
  }
  invisible(x)
}

#' @export
print.phylogp_component <- function(x, ...) print(as_kernel(x), ...)

# --- helpers for bounded kernel parameters ---------------------------------

bounded <- function(lower, upper) list(lower = lower, upper = upper)

resolve_bounds <- function(b, ti, cmp) {
  lo <- if (is.function(b$lower)) b$lower(ti, cmp) else b$lower
  hi <- if (is.function(b$upper)) b$upper(ti, cmp) else b$upper
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi)
    stop("Invalid parameter bounds.", call. = FALSE)
  c(lo, hi)
}

# Wrap an expression so the component's parameters are available under their
# short names, e.g. { real alpha = ou1_alpha; <code> }
with_aliases <- function(id, cmp, code, extra = NULL) {
  al <- vapply(names(cmp$pars), function(p)
    sprintf("real %s = %s_%s;", p, id, p), "")
  paste0("{\n", paste0("    ", c(al, extra, code), collapse = "\n"), "\n  }")
}

# One-parameter branch-wise kernels (the classical tree transformations)
theta_component <- function(name, label, par, bounds, edge, r_edge,
                            default, root = NULL, r_root = NULL,
                            form = function(ti) "edge", dense = NULL,
                            r_dense = NULL, needs = character()) {
  pars <- stats::setNames(list(bounds), par)
  new_component(
    name = name, label = label, pars = pars, form = form,
    edge = edge, root = root, dense = dense, needs = needs,
    r_eval = function(ti, p, form = "edge") {
      th <- p[[par]]
      s2 <- p$sigma^2
      if (form == "edge") {
        list(w = s2 * r_edge(ti, th),
             rv = if (is.null(r_root)) 0 else s2 * r_root(ti, th))
      } else {
        list(K = s2 * r_dense(ti, th))
      }
    },
    defaults = function(ti) stats::setNames(list(1, default(ti)), c("sigma", par)),
    validity = function(ti, form) {
      if (form == "edge") {
        list(guaranteed = TRUE,
             reason = "branch-wise form with non-negative edge variances")
      } else {
        list(guaranteed = TRUE,
             reason = "exponential (completely monotone) function of patristic distance; valid on every tree")
      }
    }
  )
}

# --- classical kernels ------------------------------------------------------

#' Classical comparative-method kernels
#'
#' Each function returns one kernel component. Components are combined with
#' `+`, e.g. `k_bm() + k_white()`. Every component has its own scale
#' parameter `sigma`; its contribution to the covariance is `sigma^2` times
#' the base kernel below (`t_ij`: shared path length from the root to the
#' most recent common ancestor of tips i and j; `d_ij`: patristic distance;
#' `T`: tree height).
#'
#' * `k_bm()`: Brownian motion, `t_ij`.
#' * `k_lambda()`: Pagel's lambda. Internal branches are multiplied by
#'   `lambda` and tip branches lengthened to keep tip heights fixed.
#' * `k_kappa()`: Pagel's kappa, every branch length raised to `kappa`.
#' * `k_delta()`: Pagel's delta, node heights transformed, giving
#'   `T^(1 - delta) t_ij^delta`.
#' * `k_eb()`: early burst / ACDC, rate `exp(rate * t)`, giving
#'   `(exp(rate * t_ij) - 1) / rate`.
#' * `k_ou()`: Ornstein-Uhlenbeck with selection strength `alpha`; random
#'   (stationary) root `exp(-alpha d_ij) / (2 alpha)`, or fixed root, which
#'   multiplies this by `1 - exp(-2 alpha t_ij)`. On ultrametric trees it is
#'   fitted in branch-wise form (O(N)); on other trees through its dense
#'   closed form.
#' * `k_white()`: independent tip variance (residual, within-species
#'   variation); `k_bm() + k_white()` is Pagel's lambda on an ultrametric tree.
#'
#' @param lower,upper Bounds of the kernel parameter. They may be numbers or
#'   functions of the tree information (see [tree_info()]); the defaults match
#'   those of `phylolm`.
#' @param root For `k_ou()`, `"random"` (stationary) or `"fixed"` root.
#' @return A kernel component.
#' @examples
#' k_bm() + k_white()
#' k_ou(root = "fixed")
#' @name classical-kernels
NULL

#' @rdname classical-kernels
#' @export
k_bm <- function() {
  new_component(
    name = "bm", label = "bm()", edge = "len",
    r_eval = function(ti, p, form = "edge") list(w = p$sigma^2 * ti$len, rv = 0),
    defaults = function(ti) list(sigma = 1),
    validity = function(ti, form) list(guaranteed = TRUE,
      reason = "branch-wise form with non-negative edge variances")
  )
}

#' @rdname classical-kernels
#' @export
k_lambda <- function(lower = 1e-7, upper = 1) {
  theta_component(
    "lambda", "lambda()", "lambda", bounded(lower, upper),
    edge = "lambda * len + (1 - lambda) * (is_tip .* he)",
    r_edge = function(ti, th) th * ti$len + (1 - th) * ti$is_tip * ti$he,
    default = function(ti) 0.5)
}

#' @rdname classical-kernels
#' @export
k_kappa <- function(lower = 1e-6, upper = 1) {
  theta_component(
    "kappa", "kappa()", "kappa", bounded(lower, upper),
    edge = "pow(len, kappa)",
    r_edge = function(ti, th) ti$len^th,
    default = function(ti) 0.5)
}

#' @rdname classical-kernels
#' @export
k_delta <- function(lower = 1e-5, upper = 3) {
  theta_component(
    "delta", "delta()", "delta", bounded(lower, upper),
    edge = "(pow(he, delta) - pow(hs, delta)) * pow(Tmax, 1 - delta)",
    r_edge = function(ti, th) (ti$he^th - ti$hs^th) * ti$Tmax^(1 - th),
    default = function(ti) 1.5)
}

#' @rdname classical-kernels
#' @export
k_eb <- function(lower = function(ti, cmp) -log(1e6) / ti$Tmax, upper = 0) {
  theta_component(
    "eb", "eb()", "rate", bounded(lower, upper),
    edge = "(abs(rate) < 1e-10) ? len : (exp(rate * he) - exp(rate * hs)) / rate",
    r_edge = function(ti, th) if (abs(th) < 1e-10) ti$len else
      (exp(th * ti$he) - exp(th * ti$hs)) / th,
    default = function(ti) -3 / ti$Tmax)
}

#' @rdname classical-kernels
#' @export
k_ou <- function(root = c("random", "fixed"),
                 lower = function(ti, cmp) 1e-7 / ti$Tmax,
                 upper = function(ti, cmp) 50 / ti$Tmax) {
  root <- match.arg(root)
  ou_edge_r <- function(ti, a)
    exp(-2 * a * ti$Tmax) / (2 * a) * (exp(2 * a * ti$he) - exp(2 * a * ti$hs))
  cmp <- theta_component(
    "ou", sprintf("ou(root = \"%s\")", root), "alpha", bounded(lower, upper),
    edge = "exp(-2 * alpha * Tmax) / (2 * alpha) * (exp(2 * alpha * he) - exp(2 * alpha * hs))",
    r_edge = ou_edge_r,
    root = if (root == "random") "exp(-2 * alpha * Tmax) / (2 * alpha)",
    r_root = if (root == "random") function(ti, a) exp(-2 * a * ti$Tmax) / (2 * a),
    # the branch-wise OU form holds on ultrametric trees only
    form = function(ti) if (ti$ultrametric) "edge" else "dense",
    dense = if (root == "random") "exp(-alpha * Dtree) / (2 * alpha)" else
      "exp(-alpha * Dtree) .* (1 - exp(-2 * alpha * tC)) / (2 * alpha)",
    r_dense = if (root == "random") function(ti, a) exp(-a * tree_D(ti)) / (2 * a) else
      function(ti, a) exp(-a * tree_D(ti)) * (1 - exp(-2 * a * tree_C(ti))) / (2 * a),
    needs = if (root == "random") "Dtree" else c("Dtree", "tC"),
    default = function(ti) 2 / ti$Tmax)
  cmp$ou_root <- root
  cmp
}

#' @rdname classical-kernels
#' @param df For `k_white()`: degrees of freedom of a Student-t residual
#'   (heavy tails, robust to outlying species). `NULL` (default) gives Gaussian
#'   residuals. The t residual is a scale mixture: each species gets a latent
#'   variance multiplier with an inverse-gamma(df/2, df/2) prior, so the model
#'   needs `method = "sample"`.
#' @export
k_white <- function(df = NULL) {
  if (is.null(df)) {
    return(new_component(
      name = "white", label = "white()", white = TRUE,
      form = function(ti) "white",
      r_eval = function(ti, p, form = "white") list(t2 = rep(p$sigma^2, ti$N)),
      defaults = function(ti) list(sigma = 1),
      validity = function(ti, form) list(guaranteed = TRUE,
        reason = "independent tip variance"),
      base_var = function(ti) 1
    ))
  }
  if (!is.numeric(df) || length(df) != 1 || df <= 0) stop("`df` must be a positive number.", call. = FALSE)
  new_component(
    name = "white", label = sprintf("white(df = %g)", df), white = TRUE, latent = TRUE,
    form = function(ti) "white",
    white_expr = "ID_lam",
    stan_data = "real<lower=0> ID_df;",
    stan_pars = "vector<lower=0>[N] ID_lam;     // variance multipliers (Student-t)",
    stan_latent_prior = "ID_lam ~ inv_gamma(0.5 * ID_df, 0.5 * ID_df);",
    standata = function(ti, id, cmp) stats::setNames(list(cmp$df), paste0(id, "_df")),
    # leave-one-out: integrate the species' own variance multiplier against
    # its inverse-gamma prior by quadrature on log(lam)
    tip_loo = list(kind = "lam", setup = character(),
                   cur = "square(ID_sigma) * ID_lam[i]",
                   node = "square(ID_sigma) * lq_lam[g]"),
    draw_pars = function(p, row, id, ti) {
      p$lam <- unname(row[paste0(id, "_lam[", seq_len(ti$N), "]")])
      p
    },
    r_eval = function(ti, p, form = "white") {
      if (!is.null(p$lam)) return(list(t2 = p$sigma^2 * p$lam))
      # marginal variance of a t residual (used for plots and checks)
      list(t2 = rep(p$sigma^2 * if (df > 2) df / (df - 2) else 1, ti$N))
    },
    defaults = function(ti) list(sigma = 1),
    validity = function(ti, form) list(guaranteed = TRUE,
      reason = "independent (Student-t) tip variance"),
    base_var = function(ti) 1,
    extra = list(df = df)
  )
}
