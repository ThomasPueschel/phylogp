# Design tools: is a kernel valid on this tree, and what does it imply?

#' Check that a kernel is a valid covariance on a tree
#'
#' A kernel must be positive semi-definite on the tips of the tree. The check
#' first applies the theory (see the package vignette and Puschel, in prep.):
#'
#' * a branch-wise kernel (edge variances accumulated along the tree) is
#'   valid whenever its edge variances are non-negative;
#' * on an ultrametric tree, any non-increasing function of patristic
#'   distance is valid;
#' * on any tree, completely monotone functions of patristic distance, such
#'   as the exponential (Ornstein-Uhlenbeck) kernel, are valid;
#' * squared-exponential and Matern kernels on patristic distances are not
#'   guaranteed valid on trees with fossil tips or without a clock.
#'
#' Where validity is not guaranteed, the check evaluates the kernel over a
#' grid of parameter values (edge variances for user-defined branch-wise
#' kernels, the smallest eigenvalue for dense kernels) and reports failures.
#' A sum of valid kernels is valid.
#'
#' @inheritParams kernel_matrix
#' @param n_grid Number of parameter values tried in numerical checks.
#' @return A data frame with one row per component (class
#'   `phylogp_check`), printed as a short report.
#' @examples
#' tr <- ape::rtree(30)                       # not ultrametric
#' kernel_check(k_dist(phi = "sqexp"), tr)    # may fail
#' kernel_check(k_dist(phi = "exponential"), tr)
#' @export
kernel_check <- function(kernel, tree, n_grid = 25) {
  ti <- tree_info(tree)
  L <- kernel_layout(kernel, ti)
  rows <- lapply(seq_along(L$ids), function(j) {
    cmp <- L$kernel$components[[j]]
    v <- cmp$validity(ti, L$forms[j])
    numeric_min <- NA_real_
    n_fail <- NA_integer_
    if (!isTRUE(v$guaranteed)) {
      grid <- par_grid(cmp, ti, n_grid)
      vals <- vapply(grid, function(p) {
        P <- stats::setNames(list(p), L$ids[j])
        kp <- kernel_parts(as_kernel(cmp), ti, P)
        if (L$forms[j] == "edge") {
          min(kp$w)
        } else {
          K <- kp$K
          min(eigen(K, symmetric = TRUE, only.values = TRUE)$values) / max(diag(K))
        }
      }, numeric(1))
      numeric_min <- min(vals)
      n_fail <- sum(vals < -1e-8)
      v$guaranteed <- n_fail == 0
      v$reason <- paste0(v$reason, if (n_fail)
        sprintf("; FAILS at %d of %d parameter values tried", n_fail, length(vals)) else
        sprintf("; no failures at %d parameter values tried", length(vals)))
    }
    data.frame(component = L$ids[j], label = cmp$label, form = L$forms[j],
               valid = v$guaranteed, reason = v$reason,
               min_value = numeric_min, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  structure(out, class = c("phylogp_check", "data.frame"), tree = ti)
}

#' @export
print.phylogp_check <- function(x, ...) {
  ti <- attr(x, "tree")
  cat(sprintf("Kernel check on a %d-tip %s tree\n", ti$N,
              if (ti$ultrametric) "ultrametric" else "non-ultrametric"))
  for (i in seq_len(nrow(x))) {
    flag <- if (isTRUE(x$valid[i])) "valid  " else "INVALID"
    cat(sprintf("  [%s] %-10s %s (%s)\n         %s\n", flag, x$component[i],
                x$label[i], x$form[i], x$reason[i]))
  }
  if (all(x$valid)) {
    cat("The kernel is a valid covariance on this tree.\n")
  } else {
    cat("The kernel is NOT guaranteed to be a valid covariance on this tree.\n",
        "Use a branch-wise kernel or a completely monotone function of distance\n",
        "(e.g. phi = \"exponential\"); adding jitter would change the model.\n", sep = "")
  }
  invisible(x)
}

# Parameter values spread over the bounds (log scale for scale-type bounds)
par_grid <- function(cmp, ti, n) {
  base <- cmp$defaults(ti)
  if (!length(cmp$pars)) return(list(base))
  bnds <- lapply(cmp$pars, resolve_bounds, ti = ti, cmp = cmp)
  lapply(seq_len(n), function(k) {
    p <- base
    for (nm in names(bnds)) {
      b <- bnds[[nm]]
      q <- if (n == 1) 0.5 else (k - 1) / (n - 1)
      # a different ordering per parameter so the grid covers combinations
      q <- if (match(nm, names(bnds)) %% 2 == 0) (q * 7) %% 1 else q
      p[[nm]] <- if (b[1] > 0 && b[2] / b[1] > 100)
        exp(log(b[1]) + q * (log(b[2]) - log(b[1]))) else b[1] + q * (b[2] - b[1])
    }
    p
  })
}

#' Plot the correlations implied by a kernel
#'
#' Draws the tree next to the correlation matrix the kernel implies between
#' its tips, with tips in the same order, so clades appear as blocks.
#'
#' @inheritParams kernel_matrix
#' @param labels Show tip labels? By default only for trees with at most 40
#'   tips.
#' @param main Title.
#' @param show_tree Draw the tree next to the matrix? With `FALSE` only the
#'   matrix is drawn, into the current plotting region (useful for panels).
#' @return The correlation matrix, invisibly.
#' @examples
#' tr <- ape::rcoal(24)
#' kernel_plot(k_bm() + k_white(), tr)
#' kernel_plot(k_ou(), tr, pars = list(ou1 = list(alpha = 5)))
#' @export
kernel_plot <- function(kernel, tree, pars = NULL, labels = NULL, main = NULL,
                        show_tree = TRUE) {
  ti <- tree_info(tree)
  K <- kernel_matrix(kernel, ti, pars)
  R <- stats::cov2cor(K)
  tr <- ape::reorder.phylo(ti$tree, "cladewise")
  ord <- tr$tip.label[tr$edge[tr$edge[, 2] <= ti$N, 2]]
  R <- R[ord, ord]
  labels <- labels %||% (ti$N <= 40)
  if (show_tree) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op))
    graphics::layout(matrix(1:2, 1), widths = c(1, 1.6))
    graphics::par(mar = c(1, 1, 3, 0))
    ape::plot.phylo(tr, show.tip.label = labels, cex = 0.6, no.margin = FALSE,
                    y.lim = c(0.5, ti$N + 0.5))
    graphics::par(mar = c(1, 0, 3, 4))
  }
  cols <- grDevices::hcl.colors(50, "YlGnBu", rev = TRUE)
  # tip k of `ord` is drawn at height k by plot.phylo, so row/column k of R
  # lines up with it
  graphics::image(seq_len(ti$N), seq_len(ti$N), R,
                  zlim = c(min(0, min(R)), 1), col = cols, axes = FALSE,
                  xlab = "", ylab = "", useRaster = TRUE)
  graphics::box()
  graphics::title(main %||% paste(vapply(as_kernel(kernel)$components, `[[`, "", "label"),
                                  collapse = " + "), cex.main = 0.95, font.main = 1,
                  line = if (show_tree) 1.5 else 0.6)
  if (show_tree)
    graphics::mtext("correlation between tips (light 0, dark 1)", side = 3,
                    line = 0.3, cex = 0.7)
  invisible(R)
}
