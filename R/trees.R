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
#'   `loo()` gives leave-one-out cross-validation with trees weighted equally;
#'   `coef()` gives per-tree estimates.
#' @references
#' Villemereuil P. de, Wells J.A., Edwards R.D., Blomberg S.P. 2012. Bayesian models for comparative analysis integrating phylogenetic uncertainty. BMC Evolutionary Biology. 12:102.
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
  # one set of regimes for every stochastic map, so the fits share parameters
  st <- sort(unique(unlist(lapply(trees[idx], function(t) colnames(t$mapped.edge)))))
  if (length(st)) kernel <- resolve_regimes(as_kernel(kernel), states = st)
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
  # PSIS-LOO on each tree separately. Running PSIS once on the draws of all
  # trees would weight each tree by the inverse of its predictive density for
  # the left-out species (a harmonic mean across trees), not equally.
  per_tree <- lapply(x$fits, function(f) {
    m <- unclass(f$fit$draws("log_lik", format = "draws_array"))
    m <- m[, , match(common, f$tree$tip.label), drop = FALSE]
    loo::loo(m, r_eff = loo::relative_eff(exp(m)), ...)
  })
  n <- length(common)
  pw <- function(col) matrix(vapply(per_tree, function(l) unname(l$pointwise[, col]), numeric(n)), n)
  elpd_t <- pw("elpd_loo")
  lpd_t <- elpd_t + pw("p_loo")
  mcse_t <- pw("mcse_elpd_loo")
  k_t <- pw("influence_pareto_k")
  # equal weights: the left-out density is the mean over trees of the
  # per-tree left-out densities
  lme <- function(a) apply(a, 1, function(v) max(v) + log(mean(exp(v - max(v)))))
  elpd <- lme(elpd_t)
  lpd <- lme(lpd_t)
  w <- exp(elpd_t - elpd) / ncol(elpd_t)
  mcse <- sqrt(rowSums((w * mcse_t)^2))
  pareto_k <- apply(k_t, 1, max)
  pointwise <- cbind(elpd_loo = elpd, mcse_elpd_loo = mcse, p_loo = lpd - elpd,
                     looic = -2 * elpd, influence_pareto_k = pareto_k)
  est <- function(v) c(sum(v), sqrt(n * stats::var(v)))
  estimates <- rbind(elpd_loo = est(elpd), p_loo = est(lpd - elpd), looic = est(-2 * elpd))
  colnames(estimates) <- c("Estimate", "SE")
  out <- per_tree[[1]]
  out$estimates <- estimates
  out$pointwise <- pointwise
  out$diagnostics$pareto_k <- pareto_k
  diag <- function(nm) matrix(vapply(per_tree, function(l) l$diagnostics[[nm]], numeric(n)), n)
  out$diagnostics$n_eff <- apply(diag("n_eff"), 1, min)
  out$diagnostics$r_eff <- rowMeans(diag("r_eff"))
  out$psis_object <- NULL
  out$elpd_loo <- estimates["elpd_loo", 1]
  out$p_loo <- estimates["p_loo", 1]
  out$looic <- estimates["looic", 1]
  out$se_elpd_loo <- estimates["elpd_loo", 2]
  out$se_p_loo <- estimates["p_loo", 2]
  out$se_looic <- estimates["looic", 2]
  attr(out, "dims") <- c(sum(vapply(per_tree, function(l) attr(l, "dims")[1], 0)), n)
  attr(out, "species") <- common
  attr(out, "per_tree") <- per_tree
  out
}
