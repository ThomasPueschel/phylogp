# Products with a covariate: phylogenetically varying regression slopes.

#' Phylogenetically varying slopes
#'
#' `k_slope(x) * k` multiplies a kernel `k` by `x_i x_j`, the kernel of a
#' regression slope on `x` whose deviations from the mean slope evolve along
#' the tree under `k`:
#'
#' \deqn{y_i = \beta_0 + (\beta_1 + v_i) x_i + \dots, \quad v \sim GP(0, k).}
#'
#' * `k_slope(x) * k_bm()` (or another branch-wise kernel): the slope evolves
#'   along the tree. The covariance is `sigma^2 (x x') o K`, valid whenever `K`
#'   is (Schur product theorem).
#' * `k_slope(x) * k_white()`: independent species-specific slope deviations,
#'   i.e. residual variance that grows with `x^2`. This is the natural control
#'   for a varying-slope model: a varying slope also makes the variance of `y`
#'   grow with distance from the centre of `x`, and comparing the two kernels
#'   separates evolving slopes from heteroscedasticity.
#'
#' The mean slope is the regression coefficient of `x` in the model formula,
#' so `x` should normally appear there as well. Because the kernel depends on
#' where `x` is centred, `x` is centred at its mean by default; report the
#' centring with the results.
#'
#' With `latent = TRUE` (default) the slope deviations at the nodes of the tree
#' are sampled as latent variables and the remaining Gaussian components are
#' integrated out by pruning, so the cost per gradient stays linear in the
#' number of species; this needs `method = "sample"`. With `latent = FALSE`
#' the dense marginal kernel is used (O(N^3)), which also allows maximum
#' likelihood.
#'
#' @param x Name of a column of the data passed to [phylo_gp()], or a numeric
#'   vector named by tip label.
#' @param center Centre `x` at its mean over the fitted species?
#' @param latent Sample the slope deviations at the nodes (linear time) rather
#'   than use the dense marginal kernel.
#' @return An object to be multiplied by one kernel component.
#' @examples
#' k_bm() + k_slope("body") * k_bm() + k_white()      # varying slopes
#' k_bm() + k_slope("body") * k_white() + k_white()   # heteroscedastic control
#' @references
#' Lande R. 1979. Quantitative genetic analysis of multivariate evolution, applied to brain:body size allometry. Evolution. 33:402–416.
#'
#' Pagel M.D., Harvey P.H. 1989. Taxonomic differences in the scaling of brain on body weight among mammals. Science. 244:1589–1593.
#'
#' Hadfield J.D., Nakagawa S. 2010. General quantitative genetic methods for comparative biology: Phylogenies, taxonomies and multi-trait models for continuous and categorical characters. Journal of Evolutionary Biology. 23:494–508.
#'
#' Matsuba M., Fukasawa K., Aoki S., Akasaka M., Ishihama F. 2024. Scalable phylogenetic Gaussian process models improve the detectability of environmental signals on local extinctions for many Red List species. Methods in Ecology and Evolution. 15:756–768.
#' @export
k_slope <- function(x, center = TRUE, latent = TRUE) {
  if (!(is.character(x) && length(x) == 1) && !(is.numeric(x) && !is.null(names(x))))
    stop("`x` must be a column name or a numeric vector named by tip label.", call. = FALSE)
  structure(list(x = x, center = center, latent = latent,
                 x_label = if (is.character(x)) x else "x"),
            class = "phylogp_scale")
}

#' @export
print.phylogp_scale <- function(x, ...) {
  cat(sprintf("k_slope(%s): multiply by a kernel component, e.g. k_slope(\"%s\") * k_bm()\n",
              x$x_label, x$x_label))
  invisible(x)
}

kernel_product <- function(e1, e2) {
  if (inherits(e1, "phylogp_scale") && !inherits(e2, "phylogp_scale")) {
    sc <- e1; base <- e2
  } else if (inherits(e2, "phylogp_scale") && !inherits(e1, "phylogp_scale")) {
    sc <- e2; base <- e1
  } else {
    stop("Products are defined between k_slope() and one kernel component.", call. = FALSE)
  }
  if (inherits(base, "phylogp_kernel")) {
    if (length(base$components) != 1)
      stop("k_slope() multiplies a single component; use (k_slope(x) * k_bm()) + ...", call. = FALSE)
    base <- base$components[[1]]
  }
  if (!inherits(base, "phylogp_component")) stop("Not a kernel component.", call. = FALSE)
  slope_component(sc, base)
}

#' Multiply a kernel by a covariate
#'
#' See [k_slope()].
#' @param e1,e2 A [k_slope()] term and one kernel component.
#' @return A kernel component.
#' @references
#' Rasmussen C.E., Williams C.K.I. 2006. Gaussian processes for machine learning. Cambridge, MA: MIT Press.
#'
#' Duvenaud D., Lloyd J., Grosse R., Tenenbaum J., Ghahramani Z. 2013. Structure discovery in nonparametric regression through compositional kernel search. Proceedings of the 30th international conference on machine learning. 28:1166–1174.
#'
#' Hadfield J.D., Nakagawa S. 2010. General quantitative genetic methods for comparative biology: Phylogenies, taxonomies and multi-trait models for continuous and categorical characters. Journal of Evolutionary Biology. 23:494–508.
#' @name kernel-product
#' @export
"*.phylogp_scale" <- kernel_product

#' @rdname kernel-product
#' @export
"*.phylogp_component" <- kernel_product

#' @rdname kernel-product
#' @export
"*.phylogp_kernel" <- kernel_product

# Values of the covariate in tip order, centred
slope_values <- function(cmp, ti) {
  x <- cmp$scale$x
  if (is.character(x))
    stop(sprintf("k_slope(\"%s\") refers to a data column; fit with phylo_gp(), or give ",
                 "a numeric vector named by tip label.", x), call. = FALSE)
  miss <- setdiff(ti$labels, names(x))
  if (length(miss)) stop("k_slope(): no value of x for ", miss[1], call. = FALSE)
  v <- as.numeric(x[ti$labels])
  if (anyNA(v)) stop("k_slope(): missing values of x.", call. = FALSE)
  centre <- cmp$scale$center_value %||% (if (cmp$scale$center) mean(v) else 0)
  v - centre
}

slope_component <- function(sc, base) {
  lab <- sprintf("slope(%s) * %s", sc$x_label, base$label)
  common <- list(stan_data = "vector[N] ID_xs;",
                 standata = function(ti, id, cmp)
                   stats::setNames(list(slope_values(cmp, ti)), paste0(id, "_xs")),
                 extra = list(scale = sc, base = base))
  if (identical(base$name, "white")) {
    if (!is.null(base$df)) stop("k_slope() * k_white(df = ) is not supported.", call. = FALSE)
    return(do.call(new_component, c(list(
      name = "slopewhite", label = lab, form = function(ti) "white",
      white_expr = "square(ID_xs)",
      r_eval = function(ti, p, form = "white") list(t2 = p$sigma^2 * slope_values(p$cmp, ti)^2),
      defaults = function(ti) list(sigma = 1),
      validity = function(ti, form) list(guaranteed = TRUE,
        reason = "independent tip variance proportional to x^2"),
      base_var = function(ti) 1), common)))
  }
  if (base$latent || is.null(base$edge) || base$name %in% c("ou", "dist"))
    stop("k_slope() can multiply k_white() or a branch-wise kernel without latent rates ",
         "(k_bm(), k_lambda(), k_kappa(), k_delta(), k_eb(), k_edge()).", call. = FALSE)
  marginal <- function(ti, p, form) {
    ev <- base$r_eval(ti, p, form = "edge")
    xs <- slope_values(p$cmp, ti)
    list(K = outer(xs, xs) * edge_matrix(ti, ev$w, ev$rv %||% 0, rep(0, ti$N)))
  }
  validity <- function(ti, form) {
    v <- base$validity(ti, "edge")
    if (isTRUE(v$guaranteed))
      v$reason <- "product of a valid branch-wise kernel and x x' (Schur product theorem)"
    v
  }
  alias <- vapply(names(base$pars), function(p) sprintf("  real %s = ID_%s;", p, p), "")
  if (!sc$latent) {
    root <- base$root %||% "0"
    return(do.call(new_component, c(list(
      name = paste0("slope", base$name), label = lab, pars = base$pars,
      form = function(ti) "dense",
      dense = sprintf("(ID_xs * ID_xs') .* tree_kernel(B, %s, %s, rep_vector(0, N))",
                      base$edge, root),
      r_eval = marginal, defaults = base$defaults, validity = validity), common)))
  }
  do.call(new_component, c(list(
    name = paste0("slope", base$name), label = lab, pars = base$pars,
    form = function(ti) "slope", latent = TRUE,
    stan_pars = "vector[E] ID_eta;",
    stan_tpars = c(
      "vector[E] ID_wb;                 // edge variances of the slope process",
      "vector[n_nodes] ID_v;            // slope deviation at each node",
      "{", alias, sprintf("  ID_wb = %s;", base$edge), "}",
      "ID_v[root] = 0;",
      "for (k in 1:E) {",
      "  int e = E - k + 1;",
      "  ID_v[des[e]] = ID_v[anc[e]] + ID_sigma * sqrt(ID_wb[e]) * ID_eta[e];",
      "}"),
    stan_latent_prior = "ID_eta ~ std_normal();",
    stan_shift = "r -= ID_xs .* ID_v[1:N];",
    tip_loo = list(kind = "closed",
                   r = "ID_xs[i] * (ID_v[i] - ID_v[anc[e]])",
                   v = "square(ID_xs[i] * ID_sigma) * ID_wb[e]"),
    draw_pars = function(p, row, id, ti) {
      p$v <- unname(row[paste0(id, "_v[", seq_len(ti$n_nodes), "]")])
      p
    },
    r_eval = function(ti, p, form = "slope") {
      if (!is.null(p$v)) return(list(shift = slope_values(p$cmp, ti) * p$v[seq_len(ti$N)]))
      marginal(ti, p, form)
    },
    defaults = base$defaults, validity = validity), common))
}

# Replace column names in k_slope() terms by the data values (named by
# species) and fix the centring at the fitted species. With several rows per
# species (`rows`), the species value is the mean over its rows, the centre is
# the mean over rows, and the value of each row is kept in `x_obs`.
resolve_slopes <- function(kernel, data, rows = NULL, species = NULL) {
  kernel$components <- lapply(kernel$components, function(cmp) {
    if (is.null(cmp$scale)) return(cmp)
    sc <- cmp$scale
    if (!is.null(rows)) {
      if (!is.character(sc$x) || !sc$x %in% names(rows))
        stop("With several rows per species, k_slope() needs a column of the data.", call. = FALSE)
      sc$x_obs <- as.numeric(rows[[sc$x]])
      m <- tapply(sc$x_obs, as.character(rows[[species]]), mean)
      sc$x <- stats::setNames(as.numeric(m), names(m))
      sc$center_value <- if (sc$center) mean(sc$x_obs) else 0
      cmp$scale <- sc
      return(cmp)
    }
    if (is.character(sc$x)) {
      if (!sc$x %in% names(data)) stop("k_slope(): no column ", sc$x, " in the data.", call. = FALSE)
      sc$x <- stats::setNames(data[[sc$x]], rownames(data))
    }
    vals <- sc$x[rownames(data)]
    sc$center_value <- if (sc$center) mean(vals) else 0
    cmp$scale <- sc
    cmp
  })
  kernel
}

# Add covariate values of new species (for prediction with dense slope kernels)
extend_slopes <- function(kernel, newdata, targets) {
  kernel$components <- lapply(kernel$components, function(cmp) {
    if (is.null(cmp$scale)) return(cmp)
    xn <- cmp$scale$x_label
    if (xn %in% names(newdata)) {
      new <- stats::setNames(newdata[[xn]], targets)
      cmp$scale$x <- c(cmp$scale$x[setdiff(names(cmp$scale$x), targets)], new)
    }
    cmp
  })
  kernel
}

slope_columns <- function(kernel) {
  unlist(lapply(kernel$components, function(cmp)
    if (!is.null(cmp$scale) && is.character(cmp$scale$x)) cmp$scale$x))
}
