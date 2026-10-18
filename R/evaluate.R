# Evaluate kernels in R: covariance matrices, simulation, likelihoods.

# Parameter values for every component: defaults, overridden by `pars`
# (a list named by component id, e.g. list(lambda1 = list(lambda = 0.3))).
resolve_pars <- function(L, ti, pars = NULL) {
  out <- vector("list", length(L$ids))
  names(out) <- L$ids
  unknown <- setdiff(names(pars), L$ids)
  if (length(unknown))
    stop("Unknown component id in `pars`: ", unknown[1],
         ". Ids are: ", paste(L$ids, collapse = ", "), call. = FALSE)
  for (j in seq_along(L$ids)) {
    cmp <- L$kernel$components[[j]]
    p <- utils::modifyList(cmp$defaults(ti), pars[[L$ids[j]]] %||% list())
    if (identical(cmp$name, "dist") && (is.null(p$ell) || is.na(p$ell))) {
      Dm <- cmp_distance(ti, cmp)
      p$ell <- stats::median(Dm[upper.tri(Dm)])
    }
    p$cmp <- cmp
    out[[j]] <- p
  }
  out
}

# Edge variances, root variance, tip variances and dense part of a kernel
kernel_parts <- function(kernel, ti, pars = NULL) {
  L <- kernel_layout(kernel, ti)
  P <- resolve_pars(L, ti, pars)
  w <- numeric(ti$E); rv <- 0; t2 <- numeric(ti$N); Kd <- NULL
  for (j in seq_along(L$ids)) {
    cmp <- L$kernel$components[[j]]
    ev <- cmp$r_eval(ti, P[[j]], form = L$forms[j])
    if (!is.null(ev$w)) w <- w + ev$w
    if (!is.null(ev$rv)) rv <- rv + ev$rv
    if (!is.null(ev$t2)) t2 <- t2 + ev$t2
    if (!is.null(ev$K)) Kd <- if (is.null(Kd)) ev$K else Kd + ev$K
  }
  list(w = w, rv = rv, t2 = t2, K = Kd, layout = L, pars = P)
}

#' Covariance matrix implied by a kernel
#'
#' @param kernel A kernel, e.g. `k_bm() + k_white()`.
#' @param tree A `phylo` tree (or the output of [tree_info()]).
#' @param pars Optional parameter values, as a list named by component id
#'   (see `print(kernel)`), e.g. `list(lambda1 = list(sigma = 1, lambda = 0.3))`.
#'   Unspecified values take the defaults of each component. For components
#'   with latent rates (`k_relaxed()`, `k_rate_time()`), a rate path is drawn
#'   from the prior unless supplied (`z` or `lograte`).
#' @return An N x N covariance matrix with tip labels as dimnames.
#' @examples
#' tr <- ape::rcoal(10)
#' kernel_matrix(k_bm() + k_white(), tr)
#' @export
kernel_matrix <- function(kernel, tree, pars = NULL) {
  ti <- tree_info(tree)
  kp <- kernel_parts(kernel, ti, pars)
  K <- edge_matrix(ti, kp$w, kp$rv, kp$t2)
  if (!is.null(kp$K)) K <- K + kp$K
  dimnames(K) <- list(ti$labels, ti$labels)
  K
}

# K = B diag(w) B' + rv + diag(t2)
edge_matrix <- function(ti, w, rv, t2) {
  B <- tree_B(ti)
  K <- B %*% (w * t(B)) + rv
  diag(K) <- diag(K) + t2
  K
}

#' Simulate traits from a kernel
#'
#' @inheritParams kernel_matrix
#' @param nsim Number of simulated data sets.
#' @param mean Mean of the trait (a number or a vector in tip order).
#' @return An N x nsim matrix with tip labels as row names.
#' @examples
#' tr <- ape::rcoal(20)
#' y <- kernel_simulate(k_ou(), tr, pars = list(ou1 = list(alpha = 3)))
#' @export
kernel_simulate <- function(kernel, tree, pars = NULL, nsim = 1, mean = 0) {
  K <- kernel_matrix(kernel, tree, pars)
  L <- chol(K + diag(1e-10 * mean(diag(K)), nrow(K)))
  Y <- mean + t(L) %*% matrix(stats::rnorm(nrow(K) * nsim), nrow(K), nsim)
  rownames(Y) <- rownames(K)
  Y
}

#' Log-likelihood of a kernel in R
#'
#' Evaluates the Gaussian log-likelihood at given parameter values, with the
#' O(N) pruning recursion when the kernel is branch-wise and a Cholesky
#' factorisation otherwise. Useful for checking a kernel or a fit.
#'
#' @inheritParams kernel_matrix
#' @param y Trait values, named by tip label or in tip order.
#' @param mean Mean vector (or number) to subtract.
#' @return The log-likelihood.
#' @export
kernel_loglik <- function(kernel, tree, y, pars = NULL, mean = 0) {
  ti <- tree_info(tree)
  r <- order_by_tips(y, ti) - mean
  kp <- kernel_parts(kernel, ti, pars)
  if (is.null(kp$K)) return(prune_loglik(ti, r, kp$w, kp$rv, kp$t2))
  K <- edge_matrix(ti, kp$w, kp$rv, kp$t2) + kp$K
  U <- chol(K)
  z <- backsolve(U, r, transpose = TRUE)
  -0.5 * sum(z^2) - sum(log(diag(U))) - 0.5 * length(r) * log(2 * pi)
}

order_by_tips <- function(y, ti) {
  if (!is.null(names(y))) {
    miss <- setdiff(ti$labels, names(y))
    if (length(miss)) stop("`y` lacks values for ", length(miss), " tips.", call. = FALSE)
    return(as.numeric(y[ti$labels]))
  }
  if (length(y) != ti$N) stop("`y` must have one value per tip.", call. = FALSE)
  as.numeric(y)
}

# O(N) pruning log-likelihood, mirrors tree_pruning_lpdf() in Stan
prune_loglik <- function(ti, r, w, rv, t2) {
  m <- numeric(ti$n_nodes); v <- numeric(ti$n_nodes)
  m[seq_len(ti$N)] <- r; v[seq_len(ti$N)] <- t2
  ll <- 0
  for (e in seq_len(ti$E)) {
    c <- ti$des[e]; p <- ti$anc[e]
    vc <- v[c] + w[e]
    if (ti$first_child[e] == 1) {
      m[p] <- m[c]; v[p] <- vc
    } else {
      s <- v[p] + vc
      ll <- ll + stats::dnorm(m[p] - m[c], 0, sqrt(s), log = TRUE)
      m[p] <- (m[p] * vc + m[c] * v[p]) / s
      v[p] <- v[p] * vc / s
    }
  }
  ll + stats::dnorm(0, m[ti$root], sqrt(v[ti$root] + rv), log = TRUE)
}

# O(N) leave-one-out predictive log densities, mirrors tree_loo_moments()
prune_loo_lpd <- function(ti, r, w, rv, t2) {
  N <- ti$N; nn <- ti$n_nodes
  um <- numeric(nn); uv <- numeric(nn); P <- numeric(nn); H <- numeric(nn)
  um[1:N] <- r; uv[1:N] <- t2
  for (e in seq_len(ti$E)) {
    c <- ti$des[e]; p <- ti$anc[e]
    if (c > N) { uv[c] <- 1 / P[c]; um[c] <- H[c] / P[c] }
    vv <- uv[c] + w[e]
    P[p] <- P[p] + 1 / vv; H[p] <- H[p] + um[c] / vv
  }
  om <- numeric(nn); ov <- numeric(nn)
  ov[ti$root] <- rv
  for (e in rev(seq_len(ti$E))) {
    c <- ti$des[e]; p <- ti$anc[e]
    vv <- uv[c] + w[e]
    if (ov[p] <= 0) {
      m_ex <- om[p]; v_ex <- 0
    } else {
      Pp <- P[p] - 1 / vv + 1 / ov[p]
      Hp <- H[p] - um[c] / vv + om[p] / ov[p]
      m_ex <- Hp / Pp; v_ex <- 1 / Pp
    }
    om[c] <- m_ex; ov[c] <- v_ex + w[e]
  }
  stats::dnorm(r, om[1:N], sqrt(ov[1:N] + t2), log = TRUE)
}
