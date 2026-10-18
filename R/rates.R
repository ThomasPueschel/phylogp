# Summaries of the latent rates estimated by k_relaxed() and k_rate_time().

find_component <- function(fit, name) {
  if (!inherits(fit, "phylogp_fit")) stop("`fit` must be a phylogp_fit.", call. = FALSE)
  if (fit$method != "sample") stop("Rate summaries need a Bayesian fit.", call. = FALSE)
  nm <- vapply(fit$kernel$components, `[[`, "", "name")
  j <- which(nm == name)
  if (!length(j))
    stop(sprintf("The kernel has no %s component.", name), call. = FALSE)
  fit$layout$ids[j[1]]
}

#' Branch rates from a relaxed-rate fit
#'
#' Posterior summaries of the rate of evolution on every edge, from a model
#' with a [k_relaxed()] component. Rates are relative to the median edge rate
#' of each posterior draw, so 2 means twice as fast as a typical branch.
#'
#' @param fit A Bayesian `phylogp_fit` whose kernel contains `k_relaxed()`.
#' @param probs Quantiles to report.
#' @return A data frame with one row per edge, in the edge order of
#'   `fit$tree`, giving the parent and child node, the posterior median and
#'   the requested quantiles of the relative rate.
#' @export
edge_rates <- function(fit, probs = c(0.05, 0.95)) {
  id <- find_component(fit, "relaxed")
  ti <- tree_info(fit$tree)
  z <- posterior::as_draws_matrix(fit$fit$draws(paste0(id, "_z")))
  lr <- 0.5 * (z[, ti$anc, drop = FALSE] + z[, ti$des, drop = FALSE])
  rel <- exp(lr - apply(lr, 1, stats::median))
  q <- t(apply(rel, 2, stats::quantile, probs = c(0.5, probs)))
  out <- data.frame(edge = seq_len(ti$E), parent = ti$anc, child = ti$des,
                    median = q[, 1])
  for (k in seq_along(probs)) out[[sprintf("q%g", 100 * probs[k])]] <- q[, k + 1]
  out
}

#' Plot branch rates on the tree
#'
#' Colours every branch by its posterior median relative rate (log scale) from
#' a relaxed-rate fit: blue slower, red faster than a typical branch.
#'
#' @param fit A Bayesian `phylogp_fit` whose kernel contains `k_relaxed()`.
#' @param ... Further arguments to [plot_tree_rates()].
#' @return The output of [edge_rates()], invisibly.
#' @export
plot_edge_rates <- function(fit, ...) {
  er <- edge_rates(fit)
  plot_tree_rates(fit$tree, er$median, ...)
  invisible(er)
}

#' Plot a tree with branches coloured by rate
#'
#' @param tree A `phylo` tree.
#' @param rates Positive rates, one per edge in the order of `tree$edge`.
#' @param lim Colour range on the log2 scale (rates from `2^-lim` to `2^lim`);
#'   by default symmetric around 1 and wide enough for all rates.
#' @param main Title.
#' @param ... Further arguments to [plot_tree_values()].
#' @return `NULL`, invisibly.
#' @export
plot_tree_rates <- function(tree, rates, lim = NULL,
                            main = "Relative rate of evolution", ...) {
  lr <- log2(rates)
  lim <- lim %||% max(1, ceiling(max(abs(lr))))
  plot_tree_values(tree, lr, center = 0, lim = lim, main = main,
                   key_labels = c(paste0("1/", signif(2^lim, 2), "x"), "1x",
                                  paste0(signif(2^lim, 2), "x")), ...)
}

#' Plot a tree with branches coloured by a value
#'
#' Colours every branch on a diverging scale (blue below `center`, red
#' above), e.g. by estimated slopes from [varying_slopes()] or by ancestral
#' states.
#'
#' @param tree A `phylo` tree.
#' @param values One value per edge, in the order of `tree$edge`.
#' @param center Value shown in the neutral colour.
#' @param lim Half-width of the colour range around `center`; values beyond
#'   it get the extreme colours. Default: the largest deviation.
#' @param type Tree layout, passed to `ape::plot.phylo()`.
#' @param show.tip.label Show tip labels? Default: only for small trees.
#' @param edge.width Branch width.
#' @param main Title.
#' @param legend Draw a colour key?
#' @param key_labels Labels for the low end, centre and high end of the key.
#' @param ... Further arguments to `ape::plot.phylo()`.
#' @return `NULL`, invisibly.
#' @export
plot_tree_values <- function(tree, values, center = stats::median(values), lim = NULL,
                             type = "phylogram", show.tip.label = NULL,
                             edge.width = 2.5, main = NULL, legend = TRUE,
                             key_labels = NULL, ...) {
  if (length(values) != nrow(tree$edge)) stop("`values` needs one value per edge.", call. = FALSE)
  dev <- values - center
  lim <- lim %||% max(abs(dev))
  if (lim <= 0) lim <- 1
  pal <- grDevices::colorRampPalette(c("#24527A", "#7FA7C9", "#ECE7DF",
                                       "#E39A7A", "#A8322A"))(101)
  col <- pal[round(50 * pmax(-1, pmin(1, dev / lim))) + 51]
  # per-edge widths and line types: ape needs them as long as the colours
  ape::plot.phylo(tree, type = type, edge.color = col,
                  edge.width = rep_len(edge.width, nrow(tree$edge)),
                  edge.lty = rep_len(1, nrow(tree$edge)),
                  show.tip.label = show.tip.label %||% (ape::Ntip(tree) <= 40),
                  cex = 0.6, main = main, cex.main = 1, font.main = 1, ...)
  if (legend) {
    usr <- graphics::par("usr")
    x0 <- usr[1] + 0.03 * diff(usr[1:2]); x1 <- usr[1] + 0.33 * diff(usr[1:2])
    y0 <- usr[3] + 0.05 * diff(usr[3:4])
    xs <- seq(x0, x1, length.out = 102)
    graphics::rect(xs[-102], y0, xs[-1], y0 + 0.025 * diff(usr[3:4]), col = pal, border = NA)
    labs <- key_labels %||% signif(c(center - lim, center, center + lim), 2)
    graphics::text(c(x0, (x0 + x1) / 2, x1), y0, labels = labs, pos = 1, cex = 0.7, xpd = TRUE)
  }
  invisible(NULL)
}

#' Species- and node-level slopes from a varying-slope fit
#'
#' Posterior summaries of the regression slope at every node of the tree,
#' `beta + v`, from a model with a latent `k_slope(x) * k` term. The mean slope
#' `beta` is the coefficient of `x` in the model formula.
#'
#' @param fit A Bayesian `phylogp_fit` with a `k_slope(x) * ...` term fitted
#'   with `latent = TRUE`.
#' @param probs Quantiles to report.
#' @return A data frame with one row per node of `fit$tree` (tips first):
#'   node number, tip label (`NA` for internal nodes), posterior mean, standard
#'   deviation and quantiles of the slope.
#' @export
varying_slopes <- function(fit, probs = c(0.05, 0.95)) {
  if (!inherits(fit, "phylogp_fit") || fit$method != "sample")
    stop("`fit` must be a Bayesian phylogp_fit.", call. = FALSE)
  j <- which(fit$layout$forms == "slope")
  if (!length(j)) stop("The kernel has no latent k_slope() term.", call. = FALSE)
  j <- j[1]
  id <- fit$layout$ids[j]
  xname <- fit$kernel$components[[j]]$scale$x_label
  k <- match(xname, fit$coef_names)
  if (is.na(k)) stop("`", xname, "` is not in the model formula, so the mean slope is 0; ",
                     "add it to the formula.", call. = FALSE)
  ti <- tree_info(fit$tree)
  v <- posterior::as_draws_matrix(fit$fit$draws(paste0(id, "_v")))
  b <- as.numeric(posterior::as_draws_matrix(fit$fit$draws(sprintf("coefs[%d]", k))))
  sl <- unclass(v) + b
  q <- t(apply(sl, 2, stats::quantile, probs = probs, names = FALSE))
  out <- data.frame(node = seq_len(ti$n_nodes),
                    label = c(ti$labels, rep(NA, ti$n_nodes - ti$N)),
                    mean = colMeans(sl), sd = apply(sl, 2, stats::sd))
  for (m in seq_along(probs)) out[[sprintf("q%g", 100 * probs[m])]] <- q[, m]
  rownames(out) <- NULL
  out
}

#' Rate of evolution through time
#'
#' Posterior summaries of the rate of evolution in each time bin of a
#' [k_rate_time()] component.
#'
#' @param fit A Bayesian `phylogp_fit` whose kernel contains `k_rate_time()`.
#' @param probs Quantiles to report.
#' @return A data frame with one row per bin: start and end time (from the
#'   root), posterior median and quantiles of the rate `sigma^2(t)`.
#' @export
rate_through_time <- function(fit, probs = c(0.05, 0.95)) {
  id <- find_component(fit, "ratetime")
  ti <- tree_info(fit$tree)
  M <- fit$standata[[paste0(id, "_M")]]
  f <- posterior::as_draws_matrix(fit$fit$draws(paste0(id, "_lograte")))
  s2 <- as.numeric(posterior::as_draws_matrix(fit$fit$draws(paste0(id, "_sigma"))))^2
  rate <- exp(f) * s2
  q <- t(apply(rate, 2, stats::quantile, probs = c(0.5, probs)))
  brk <- seq(0, ti$Tmax, length.out = M + 1)
  out <- data.frame(start = brk[-(M + 1)], end = brk[-1], median = q[, 1])
  for (k in seq_along(probs)) out[[sprintf("q%g", 100 * probs[k])]] <- q[, k + 1]
  out
}
