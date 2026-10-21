# Regression effects: the expected response along one predictor.

#' Plot regression effects
#'
#' Plots the expected response, \eqn{x^\top\beta}, against one predictor at a
#' time, with a posterior credible band, over the data. Other predictors are
#' held at their mean (numeric) or reference level (factors). The line is the
#' population-level relationship, averaged over the phylogenetic and other
#' kernel effects, as in a mixed model. For maximum-likelihood fits the line
#' is drawn without a band.
#'
#' Named differently from `brms::conditional_effects()` so that both
#' packages can be loaded together.
#'
#' @param fit A `phylogp_fit`.
#' @param effect Names of predictors (variables in the formula) to plot. By
#'   default every predictor, one panel each.
#' @param probs Lower and upper quantiles of the band.
#' @param points Draw the data? These are the species values, or the rows for
#'   fits with several rows per species.
#' @param n Number of grid values for a numeric predictor.
#' @param ndraws Number of posterior draws used.
#' @param plot Draw the plot? If `FALSE`, only the values are returned.
#' @param col Colour of the line and band.
#' @param ... Further graphical parameters passed to [graphics::plot()].
#' @return Invisibly, a data frame with the predictor (`effect`), its value
#'   (`value` for numeric predictors, `level` for factors), the posterior
#'   median of the expected response (`estimate`) and the band (`lower`,
#'   `upper`).
#' @examples
#' \dontrun{
#' fit <- phylo_gp(brain ~ body, d, tree, k_bm() + k_white())
#' plot_effects(fit)
#' plot_effects(fit, "body", probs = c(0.025, 0.975))
#' }
#' @export
plot_effects <- function(fit, effect = NULL, probs = c(0.05, 0.95), points = TRUE,
                         n = 100, ndraws = 1000, plot = TRUE, col = "#4F6D8F", ...) {
  if (!inherits(fit, "phylogp_fit")) stop("`fit` must be a phylogp_fit.", call. = FALSE)
  tt <- stats::delete.response(fit$terms)
  vars <- all.vars(tt)
  if (!length(vars)) stop("The model has no predictors.", call. = FALSE)
  effect <- effect %||% vars
  bad <- setdiff(effect, vars)
  if (length(bad)) stop("Not a predictor in the model: ", paste(bad, collapse = ", "), call. = FALSE)
  dat <- if (is.null(fit$repeated)) fit$data else fit$repeated$rows
  mf <- stats::model.frame(fit$terms, dat)
  y <- stats::model.response(mf)
  ylab <- deparse(fit$formula[[2]])
  D <- fit_draws(fit, ndraws)
  B <- D[, paste0("coefs[", seq_along(fit$coef_names), "]"), drop = FALSE]
  ref <- lapply(dat[vars], function(v) {
    if (is.numeric(v)) mean(v) else factor(levels(factor(v))[1], levels = levels(factor(v)))
  })

  out <- do.call(rbind, lapply(effect, function(v) {
    x <- dat[[v]]
    grid <- if (is.numeric(x)) seq(min(x), max(x), length.out = n) else levels(factor(x))
    nd <- as.data.frame(lapply(ref, function(r) rep(r, length(grid))))
    nd[[v]] <- if (is.numeric(x)) grid else factor(grid, levels = levels(factor(x)))
    X <- stats::model.matrix(tt, stats::model.frame(tt, nd, xlev = fit$xlevels), xlev = fit$xlevels)
    mu <- X %*% t(B)
    q <- t(apply(mu, 1, stats::quantile, probs = c(0.5, probs), names = FALSE))
    data.frame(effect = v,
               value = if (is.numeric(x)) grid else NA_real_,
               level = if (is.numeric(x)) NA_character_ else grid,
               estimate = q[, 1],
               lower = if (nrow(B) > 1) q[, 2] else NA_real_,
               upper = if (nrow(B) > 1) q[, 3] else NA_real_)
  }))
  rownames(out) <- NULL

  if (plot) {
    if (length(effect) > 1) {
      op <- graphics::par(mfrow = rev(grDevices::n2mfrow(length(effect))))
      on.exit(graphics::par(op))
    }
    for (v in effect) {
      e <- out[out$effect == v, ]
      x <- dat[[v]]
      band <- !all(is.na(e$lower))
      if (is.numeric(x)) {
        g <- e$value
        yl <- range(if (points) y, e$estimate, if (band) c(e$lower, e$upper), na.rm = TRUE)
        graphics::plot(range(g), yl, type = "n", las = 1, bty = "l", xlab = v, ylab = ylab, ...)
        if (points) graphics::points(x, y, pch = 19, cex = 0.7, col = grDevices::adjustcolor("grey30", 0.35))
        if (band) graphics::polygon(c(g, rev(g)), c(e$lower, rev(e$upper)),
                                    col = grDevices::adjustcolor(col, 0.3), border = NA)
        graphics::lines(g, e$estimate, lwd = 2.5, col = col)
      } else {
        lv <- e$level; k <- seq_along(lv)
        yl <- range(if (points) y, e$estimate, if (band) c(e$lower, e$upper), na.rm = TRUE)
        graphics::plot(range(k) + c(-0.5, 0.5), yl, type = "n", las = 1, bty = "l", xaxt = "n",
                       xlab = v, ylab = ylab, ...)
        graphics::axis(1, at = k, labels = lv)
        if (points) graphics::points(match(as.character(x), lv) + stats::runif(length(x), -0.15, 0.15), y,
                                     pch = 19, cex = 0.7, col = grDevices::adjustcolor("grey30", 0.35))
        if (band) graphics::segments(k, e$lower, k, e$upper, lwd = 2.5, col = col)
        graphics::points(k, e$estimate, pch = 19, cex = 1.4, col = col)
      }
    }
  }
  invisible(out)
}
