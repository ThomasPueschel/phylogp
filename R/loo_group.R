# Leave-group-out (e.g. leave-clade-out) cross-validation.

#' Leave-clade-out cross-validation
#'
#' Leave-one-out mostly asks how well a species is predicted from its closest
#' relatives. Leaving out whole clades (or any other groups, such as orders
#' or families) asks a harder question that bears more directly on the
#' evolutionary process: how well does the model predict a group of species
#' from the rest of the tree?
#'
#' For each group the model is refitted without the group's data (its
#' species are given a very large measurement error, so the tree and any
#' latent values on it stay the same), and the joint predictive density of
#' the group, `log mean_s p(y_group | y_rest, theta_s)`, is computed. For
#' branch-wise kernels this is the difference of two pruning likelihoods, so
#' it costs O(N) per draw; dense kernels use dense conditioning.
#'
#' @param fit A Bayesian `phylogp_fit`.
#' @param groups A vector giving the group of each species, named by species
#'   (tip label), or an unnamed vector in the order of the data rows used in
#'   the fit (`rownames(fit$data)`).
#' @param ndraws Number of posterior draws of each refit used.
#' @param ... Further arguments to `$sample()`.
#' @return An object of class `c("kfold", "loo")` usable with
#'   `loo::loo_compare()`: estimates of the expected log predictive density
#'   summed over groups (`elpd_kfold`), and a `pointwise` matrix with one row
#'   per group. The attribute `groups` gives group names and sizes, and
#'   `species_elpd` the marginal predictive density of each left-out species.
#' @examples
#' \dontrun{
#' lco <- loo_group(fit, groups = setNames(d$order, rownames(d)))
#' loo::loo_compare(lco_model_a, lco_model_b)
#' }
#' @export
loo_group <- function(fit, groups, ndraws = 400, ...) {
  if (!inherits(fit, "phylogp_fit") || fit$method != "sample")
    stop("`fit` must be a Bayesian phylogp_fit.", call. = FALSE)
  sp <- fit$tree$tip.label
  if (is.null(names(groups))) {
    if (length(groups) != nrow(fit$data)) stop("`groups` needs one value per species.", call. = FALSE)
    names(groups) <- rownames(fit$data)
  }
  g <- groups[sp]
  if (anyNA(g)) stop("`groups` lacks some fitted species, e.g. ", sp[is.na(g)][1], call. = FALSE)
  levels_g <- unique(as.character(g))
  ti <- tree_info(fit$tree)
  L <- kernel_layout(fit$kernel, ti)
  mod <- compile_phylogp(fit$code)
  y <- fit$standata$y; X <- fit$standata$X; se <- fit$standata$se
  ncoef <- length(fit$coef_names)
  big <- 1e4 * (stats::sd(y) + max(se))
  sm <- fit$sampler
  elpd <- numeric(length(levels_g))
  sp_elpd <- stats::setNames(numeric(length(sp)), sp)
  for (k in seq_along(levels_g)) {
    idx <- which(g == levels_g[k])
    message(sprintf("Group %d of %d: %s (%d species)", k, length(levels_g), levels_g[k], length(idx)))
    sd_k <- fit$standata
    sd_k$se[idx] <- big
    if (!is.null(fit$repeated)) sd_k$within_include[idx] <- 0
    f <- mod$sample(data = sd_k, chains = sm$chains, parallel_chains = sm$chains,
                    iter_warmup = sm$iter_warmup, iter_sampling = sm$iter_sampling,
                    adapt_delta = sm$adapt_delta, seed = sm$seed, refresh = 0,
                    show_messages = FALSE, ...)
    D <- unclass(posterior::as_draws_matrix(f$draws()))
    if (nrow(D) > ndraws) D <- D[round(seq(1, nrow(D), length.out = ndraws)), , drop = FALSE]
    res <- vapply(seq_len(nrow(D)), function(s) {
      row <- D[s, ]
      beta <- unname(row[paste0("coefs[", seq_len(ncoef), "]")])
      kp <- kernel_parts(L$kernel, ti, draw_pars(L, row, ti))
      r <- y - as.numeric(X %*% beta) - kp$shift
      t2 <- kp$t2 + se^2 + within_t2(fit, row)
      wl <- within_lp(fit, row, beta)[idx]
      r_rest <- r
      r_rest[idx] <- NA
      joint <- group_lpd(ti, kp, r, r_rest, t2, idx)
      cm <- condition_tips(ti, kp, r_rest, t2, idx)
      c(joint + sum(wl), stats::dnorm(r[idx], cm$mean, sqrt(cm$var + t2[idx]), log = TRUE) + wl)
    }, numeric(1 + length(idx)))
    res <- matrix(res, nrow = 1 + length(idx))
    lme <- function(x) max(x) + log(mean(exp(x - max(x))))
    elpd[k] <- lme(res[1, ])
    sp_elpd[idx] <- apply(res[-1, , drop = FALSE], 1, lme)
  }
  G <- length(elpd)
  pointwise <- cbind(elpd_kfold = elpd, p_kfold = NA_real_, kfoldic = -2 * elpd)
  est <- rbind(elpd_kfold = c(sum(elpd), sqrt(G * stats::var(elpd))),
               p_kfold = c(NA, NA),
               kfoldic = c(-2 * sum(elpd), 2 * sqrt(G * stats::var(elpd))))
  colnames(est) <- c("Estimate", "SE")
  out <- list(estimates = est, pointwise = pointwise)
  attr(out, "K") <- G
  attr(out, "groups") <- data.frame(group = levels_g, n = as.integer(table(g)[levels_g]),
                                    elpd = elpd)
  attr(out, "species_elpd") <- sp_elpd
  class(out) <- c("kfold", "loo")
  out
}

# Joint log density of the residuals of tips `idx` given the other tips:
# log p(r) - log p(r_rest), O(N) for branch-wise kernels.
group_lpd <- function(ti, kp, r, r_rest, t2, idx) {
  if (is.null(kp$K)) {
    return(prune_loglik_na(ti, r, kp$w, kp$rv, t2) -
             prune_loglik_na(ti, r_rest, kp$w, kp$rv, t2))
  }
  K <- edge_matrix(ti, kp$w, kp$rv, t2) + kp$K
  dmvn <- function(x, S) {
    U <- chol(S)
    z <- backsolve(U, x, transpose = TRUE)
    -0.5 * sum(z^2) - sum(log(diag(U))) - 0.5 * length(x) * log(2 * pi)
  }
  o <- which(!is.na(r_rest))
  dmvn(r, K) - dmvn(r[o], K[o, o, drop = FALSE])
}

# Pruning log-likelihood that skips unobserved tips (NA), O(N).
# Messages are kept as (mean, variance) of a Gaussian likelihood for the node
# value; a node without observed descendants carries no information.
prune_loglik_na <- function(ti, r, w, rv, t2) {
  nn <- ti$n_nodes; N <- ti$N
  m <- numeric(nn); v <- rep(Inf, nn)
  obs <- !is.na(r)
  m[seq_len(N)][obs] <- r[obs]
  v[seq_len(N)][obs] <- t2[obs]
  ll <- 0
  for (e in seq_len(ti$E)) {
    c <- ti$des[e]; p <- ti$anc[e]
    if (!is.finite(v[c])) next                     # nothing observed below c
    vc <- v[c] + w[e]
    if (!is.finite(v[p])) {
      m[p] <- m[c]; v[p] <- vc
    } else {
      s <- v[p] + vc
      ll <- ll + stats::dnorm(m[p] - m[c], 0, sqrt(s), log = TRUE)
      m[p] <- (m[p] * vc + m[c] * v[p]) / s
      v[p] <- v[p] * vc / s
    }
  }
  if (!is.finite(v[ti$root])) return(0)
  ll + stats::dnorm(0, m[ti$root], sqrt(v[ti$root] + rv), log = TRUE)
}
