# Regime-dependent rates: the rate of evolution depends on a discrete state
# (a selective regime, habitat, diet...) mapped onto the branches of the tree.

#' Rates that depend on mapped regimes
#'
#' Brownian motion whose rate depends on a discrete regime painted on the
#' tree, as in multi-rate Brownian motion models (O'Meara et al. 2006). The
#' regimes come from the tree: give [phylo_gp()] a stochastic character map
#' (a `simmap` from `phytools::make.simmap()`, or a tree painted with
#' `phytools::paintSubTree()`), whose `mapped.edge` holds the time each
#' branch spends in each regime. Regimes may change along a branch.
#'
#' Each edge accumulates variance `sigma^2 * sum_r rate_r * t_r`, where `t_r`
#' is the time it spends in regime `r`. The rate of the reference regime is
#' 1, so `sigma` is its standard deviation per unit time and the `rate_`
#' parameters are the rates of the other regimes relative to it. The kernel
#' is branch-wise, so it keeps the linear-time likelihood, leave-one-out,
#' prediction and ancestral states, and it can be fitted by maximum
#' likelihood. In Bayesian fits the relative rates get lognormal priors
#' centred on 1.
#'
#' To integrate over uncertainty in the history of the regimes, fit a sample
#' of stochastic maps (a `multiSimmap`) with [phylo_gp_trees()].
#'
#' @param states Optional regime names. By default, all regimes mapped on the
#'   tree (sorted).
#' @param reference Regime whose rate is 1; by default the first of `states`.
#' @param rate_scale Standard deviation of the lognormal priors of the
#'   relative rates.
#' @param lower,upper Bounds of the relative rates.
#' @return A kernel component.
#' @examples
#' \dontrun{
#' tr <- phytools::pbtree(n = 100, scale = 1)
#' Q <- matrix(c(-1, 1, 1, -1), 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
#' smap <- phytools::sim.history(tr, Q)
#' fit <- phylo_gp(y ~ 1, d, smap, k_regime() + k_white(), method = "optimize")
#' }
#' @references
#' O'Meara B.C., Ané C., Sanderson M.J., Wainwright P.C. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution. 60:922–933.
#'
#' Huelsenbeck J.P., Nielsen R., Bollback J.P. 2003. Stochastic mapping of morphological characters. Systematic Biology. 52:131–158.
#'
#' Revell L.J. 2012. phytools: An R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution. 3:217–223.
#' @export
k_regime <- function(states = NULL, reference = NULL, rate_scale = 1,
                     lower = 1e-3, upper = 1e3) {
  args <- list(reference = reference, rate_scale = rate_scale, lower = lower, upper = upper)
  if (!is.null(states)) return(do.call(regime_component, c(list(states = states), args)))
  # regimes are read from the tree when the kernel meets one (kernel_layout())
  new_component(name = "regime", label = "regime()",
                defaults = function(ti) list(sigma = 1),
                extra = list(regime_pending = TRUE, regime_args = args))
}

regime_component <- function(states, reference = NULL, rate_scale = 1,
                             lower = 1e-3, upper = 1e3) {
  states <- as.character(states)
  if (anyDuplicated(states)) stop("k_regime(): duplicated regime names.", call. = FALSE)
  if (length(states) < 2) stop("k_regime() needs at least two regimes.", call. = FALSE)
  reference <- as.character(reference %||% states[1])
  if (!reference %in% states) stop("k_regime(): `reference` is not a regime.", call. = FALSE)
  states <- c(reference, setdiff(states, reference))
  pn <- paste0("rate_", regime_names(states[-1]))
  edge <- paste(c("ID_T[, 1]", sprintf("%s * ID_T[, %d]", pn, seq_along(pn) + 1)),
                collapse = " + ")
  new_component(
    name = "regime", label = sprintf("regime(reference = \"%s\")", reference),
    pars = stats::setNames(rep(list(bounded(lower, upper)), length(pn)), pn),
    edge = edge,
    stan_data = c("int<lower=2> ID_R;", "matrix<lower=0>[E, ID_R] ID_T;",
                  "real<lower=0> ID_rate_scale;"),
    stan_prior = sprintf("ID_%s ~ lognormal(0, ID_rate_scale);", pn),
    standata = function(ti, id, cmp) stats::setNames(
      list(length(cmp$states), regime_time(ti, cmp$states), cmp$rate_scale),
      paste0(id, c("_R", "_T", "_rate_scale"))),
    r_eval = function(ti, p, form = "edge") {
      rates <- c(1, vapply(pn, function(nm) p[[nm]], numeric(1)))
      list(w = p$sigma^2 * as.numeric(regime_time(ti, states) %*% rates), rv = 0)
    },
    defaults = function(ti) c(list(sigma = 1), stats::setNames(as.list(rep(1, length(pn))), pn)),
    validity = function(ti, form) list(guaranteed = TRUE,
      reason = "branch-wise form; edge variances are positive rates times time in each regime"),
    extra = list(states = states, reference = reference, rate_scale = rate_scale,
                 rate_states = stats::setNames(states[-1], pn))
  )
}

# Regime names usable in Stan identifiers
regime_names <- function(s) {
  x <- gsub("[^A-Za-z0-9]+", "_", s)
  x <- gsub("^_+|_+$", "", x)
  x[x == ""] <- "r"
  make.unique(x, sep = "_")
}

# Fill in the regimes of k_regime() terms from a tree
resolve_regimes <- function(kernel, ti = NULL, states = NULL) {
  kernel$components <- lapply(kernel$components, function(cmp) {
    if (!isTRUE(cmp$regime_pending)) return(cmp)
    st <- states %||% if (!is.null(ti$regime)) sort(colnames(ti$regime))
    if (is.null(st))
      stop("k_regime() needs a tree with mapped regimes (a phytools simmap, e.g. from ",
           "make.simmap() or paintSubTree()), or its `states`.", call. = FALSE)
    do.call(regime_component, c(list(states = st), cmp$regime_args))
  })
  kernel
}

# Time each edge spends in each regime, columns in the order of `states`
regime_time <- function(ti, states) {
  M <- ti$regime
  if (is.null(M))
    stop("k_regime(): the tree has no mapped regimes (a phytools simmap is needed).", call. = FALSE)
  extra <- setdiff(colnames(M), states)
  if (length(extra) && any(M[, extra] > 0))
    stop("The tree has regimes the kernel does not know: ", paste(extra, collapse = ", "), call. = FALSE)
  out <- matrix(0, ti$E, length(states), dimnames = list(NULL, states))
  have <- intersect(states, colnames(M))
  out[, have] <- M[, have]
  out
}

# Checked mapped.edge of a tree, or NULL
mapped_edge <- function(tree) {
  M <- tree$mapped.edge
  if (is.null(M)) return(NULL)
  M <- as.matrix(M)
  if (nrow(M) != nrow(tree$edge) || is.null(colnames(M)))
    stop("`mapped.edge` must have one row per edge and regime names as column names.", call. = FALSE)
  if (any(abs(rowSums(M) - tree$edge.length) > 1e-6 * max(tree$edge.length)))
    stop("The times in `mapped.edge` do not add up to the branch lengths.", call. = FALSE)
  M
}

# ape::drop.tip() that keeps the regime map: an edge of the pruned tree may
# join several edges of the full tree, so its time in each regime is taken
# from the cumulative times (from the root) of its two end nodes.
drop_tips <- function(tree, tips) {
  M <- mapped_edge(tree)
  if (!length(tips) || is.null(M)) return(if (length(tips)) ape::drop.tip(tree, tips) else tree)
  n <- ape::Ntip(tree); nn <- n + tree$Nnode
  cum <- matrix(0, nn, ncol(M), dimnames = list(NULL, colnames(M)))
  for (e in ape::reorder.phylo(tree, "cladewise", index.only = TRUE))
    cum[tree$edge[e, 2], ] <- cum[tree$edge[e, 1], ] + M[e, ]
  labels <- tree$node.label
  tree$maps <- NULL
  tree$mapped.edge <- NULL
  tree$node.label <- paste0(".phylogp_node", n + seq_len(tree$Nnode))
  class(tree) <- "phylo"                  # not phytools' drop.tip() method for simmaps
  pr <- ape::drop.tip(tree, tips)
  orig <- c(match(pr$tip.label, tree$tip.label),
            as.integer(sub(".phylogp_node", "", pr$node.label, fixed = TRUE)))
  pr$mapped.edge <- cum[orig[pr$edge[, 2]], , drop = FALSE] - cum[orig[pr$edge[, 1]], , drop = FALSE]
  pr$mapped.edge[pr$mapped.edge < 0] <- 0
  pr$node.label <- if (is.null(labels)) NULL else labels[orig[-seq_len(ape::Ntip(pr))] - n]
  pr
}
