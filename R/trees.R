# Phylogenetic uncertainty: fitting a model across a sample of trees.

#' Fit a model across a sample of trees
#'
#' Fits the same model to each tree of a sample (for example, trees drawn
#' from the posterior of a dated phylogeny) and pools the results, so that
#' estimates and model comparisons integrate over phylogenetic uncertainty.
#' Each tree gets the same number of posterior draws, i.e. trees are weighted
#' equally, as in the common practice of repeating an analysis over a tree
#' sample. The Stan program is compiled once and reused for every tree.
#'
#' @param formula,data,kernel,species,se As in [phylo_gp()].
#' @param trees A `multiPhylo` object or a list of `phylo` trees.
#' @param n_trees Optional number of trees to use, drawn at random from
#'   `trees`; all trees by default.
#' @param method `"sample"` (pooled posterior) or `"optimize"` (the
#'   distribution of maximum-likelihood estimates across trees).
#' @param seed Random seed for tree selection and fitting.
#' @param ... Further arguments to [phylo_gp()], e.g. `chains`,
#'   `iter_sampling`.
#' @return An object of class `phylogp_trees` holding the individual fits.
#'   `summary()` gives pooled estimates and the spread of per-tree estimates;
#'   `loo()` gives leave-one-out cross-validation for the pooled posterior;
#'   `coef()` gives per-tree estimates.
#' @export
phylo_gp_trees <- function(formula, data, trees, kernel, n_trees = NULL,
                           method = c("sample", "optimize"), species = NULL,
                           se = NULL, seed = 1, ...) {
  method <- match.arg(method)
  if (inherits(trees, "phylo")) trees <- list(trees)
  trees <- unclass(trees)
  if (!length(trees) || !all(vapply(trees, inherits, NA, "phylo")))
    stop("`trees` must be a multiPhylo object or a list of phylo trees.", call. = FALSE)
  set.seed(seed)
  idx <- if (is.null(n_trees) || n_trees >= length(trees)) seq_along(trees) else
    sort(sample(length(trees), n_trees))
  fits <- vector("list", length(idx))
  for (k in seq_along(idx)) {
    message(sprintf("Tree %d of %d", k, length(idx)))
    fits[[k]] <- suppressMessages(phylo_gp(formula, data, trees[[idx[k]]], kernel,
                                           method = method, species = species, se = se,
                                           seed = seed + k, ...))
  }
  sp <- lapply(fits, function(f) sort(f$tree$tip.label))
  if (length(unique(sp)) > 1)
    warning("The trees do not all contain the same species with data; pooled ",
            "leave-one-out uses the species common to all trees.", call. = FALSE)
  structure(list(fits = fits, tree_index = idx, method = method, formula = formula,
                 kernel = fits[[1]]$kernel),
            class = "phylogp_trees")
}

#' @export
print.phylogp_trees <- function(x, digits = 3, ...) {
  cat("Phylogenetic Gaussian process fitted to", length(x$fits), "trees\n")
  cat("  Formula:", deparse(x$formula), "\n")
  cat("  Kernel: ", kernel_label(x$kernel), "\n")
  cat("  Method: ", if (x$method == "sample") "Bayesian (HMC), pooled across trees" else
    "maximum likelihood on each tree", "\n\n")
  print(summary(x), digits = digits)
  invisible(x)
}

#' @export
summary.phylogp_trees <- function(object, probs = c(0.05, 0.95), ...) {
  f1 <- object$fits[[1]]
  pn <- fit_par_names(f1)
  lab <- pn
  lab[seq_along(f1$coef_names)] <- f1$coef_names
  per_tree <- stats::coef(object)
  if (object$method == "optimize") {
    out <- data.frame(parameter = lab, mean = colMeans(per_tree),
                      between_tree_sd = apply(per_tree, 2, stats::sd),
                      min = apply(per_tree, 2, min), max = apply(per_tree, 2, max))
  } else {
    pooled <- do.call(rbind, lapply(object$fits, function(f)
      unclass(posterior::as_draws_matrix(f$fit$draws(pn)))))
    q <- t(apply(pooled, 2, stats::quantile, probs = probs, names = FALSE))
    out <- data.frame(parameter = lab, mean = colMeans(pooled),
                      sd = apply(pooled, 2, stats::sd))
    for (k in seq_along(probs)) out[[sprintf("q%g", 100 * probs[k])]] <- q[, k]
    # how much of the posterior spread comes from phylogenetic uncertainty
    out$between_tree_sd <- apply(per_tree, 2, stats::sd)
    rhat <- vapply(object$fits, function(f)
      max(f$fit$summary(pn, "rhat")$rhat, na.rm = TRUE), 0)
    attr(out, "max_rhat") <- max(rhat)
  }
  rownames(out) <- NULL
  out
}

#' @export
coef.phylogp_trees <- function(object, ...) {
  m <- t(vapply(object$fits, stats::coef, numeric(length(fit_par_names(object$fits[[1]])))))
  rownames(m) <- paste0("tree", object$tree_index)
  m
}

#' @rdname loo.phylogp_fit
#' @exportS3Method
loo.phylogp_trees <- function(x, ...) {
  if (x$method != "sample") stop("loo() needs Bayesian fits.", call. = FALSE)
  common <- Reduce(intersect, lapply(x$fits, function(f) f$tree$tip.label))
  common <- sort(common)
  ll <- lapply(x$fits, function(f) {
    m <- unclass(f$fit$draws("log_lik", format = "draws_array"))
    m[, , match(common, f$tree$tip.label), drop = FALSE]
  })
  n_iter <- min(vapply(ll, function(a) dim(a)[1], 0))
  n_chain <- vapply(ll, function(a) dim(a)[2], 0)
  # stack the trees as additional chains: iterations x chains x species
  arr <- array(NA_real_, c(n_iter, sum(n_chain), length(common)))
  at <- 0
  for (k in seq_along(ll)) {
    arr[, at + seq_len(n_chain[k]), ] <- ll[[k]][seq_len(n_iter), , , drop = FALSE]
    at <- at + n_chain[k]
  }
  # relative efficiency within each tree (draws from different trees target
  # different posteriors, so they are not chains of one sampler)
  r_eff <- Reduce(`+`, lapply(ll, function(a) loo::relative_eff(exp(a[seq_len(n_iter), , , drop = FALSE])))) / length(ll)
  out <- loo::loo(arr, r_eff = r_eff, ...)
  attr(out, "species") <- common
  out
}
