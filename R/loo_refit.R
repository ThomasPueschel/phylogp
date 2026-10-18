# Exact leave-one-out for influential species.

#' Exact leave-one-out for species with high Pareto k
#'
#' PSIS-LOO approximates leaving a species out by reweighting the posterior
#' draws of the full fit. For influential species (Pareto k above about 0.7)
#' the approximation is unreliable. `loo_refit()` refits the model without
#' each such species and replaces its contribution by the exact value
#' `log mean_s p(y_i | y_-i, theta_s)`, with `theta_s` drawn from the posterior
#' given the other species.
#'
#' A species is left out by giving it a measurement-error standard deviation
#' many times larger than the spread of the data, so its value carries no
#' information while the tree, and any latent rates on it, stay the same. Its
#' predictive density is then computed in linear time (branch-wise kernels) or
#' by dense conditioning.
#'
#' @param fit A Bayesian `phylogp_fit`.
#' @param loo Optional result of `loo(fit)`; computed if missing.
#' @param k_threshold Species with Pareto k above this value are refitted.
#' @param ... Further arguments to `$sample()`.
#' @return The `loo` object with the pointwise and total estimates updated for
#'   the refitted species (their Pareto k set to 0, as in `loo::reloo()`), and
#'   an attribute `refitted` naming them.
#' @export
loo_refit <- function(fit, loo = NULL, k_threshold = 0.7, ...) {
  if (!inherits(fit, "phylogp_fit") || fit$method != "sample")
    stop("`fit` must be a Bayesian phylogp_fit.", call. = FALSE)
  lo <- loo %||% loo.phylogp_fit(fit)
  k <- lo$diagnostics$pareto_k
  idx <- which(k > k_threshold)
  if (!length(idx)) {
    message("No species with Pareto k above ", k_threshold, ".")
    return(lo)
  }
  sp <- fit$tree$tip.label
  message(sprintf("Refitting without %d species: %s%s", length(idx),
                  paste(utils::head(sp[idx], 5), collapse = ", "),
                  if (length(idx) > 5) ", ..." else ""))
  mod <- compile_phylogp(fit$code)
  ti <- tree_info(fit$tree)
  L <- kernel_layout(fit$kernel, ti)
  y <- fit$standata$y; X <- fit$standata$X; se <- fit$standata$se
  ncoef <- length(fit$coef_names)
  big <- 1e4 * (stats::sd(y) + max(se))
  sm <- fit$sampler
  for (i in idx) {
    sd_i <- fit$standata
    sd_i$se[i] <- big
    f <- mod$sample(data = sd_i, chains = sm$chains, parallel_chains = sm$chains,
                    iter_warmup = sm$iter_warmup, iter_sampling = sm$iter_sampling,
                    adapt_delta = sm$adapt_delta, seed = sm$seed, refresh = 0,
                    show_messages = FALSE, ...)
    D <- unclass(posterior::as_draws_matrix(f$draws()))
    lp <- vapply(seq_len(nrow(D)), function(s) {
      row <- D[s, ]
      beta <- unname(row[paste0("coefs[", seq_len(ncoef), "]")])
      kp <- kernel_parts(L$kernel, ti, draw_pars(L, row, ti))
      r <- y - as.numeric(X %*% beta)
      t2 <- kp$t2 + se^2
      r_i <- r[i]
      r[i] <- NA
      cm <- condition_tips(ti, kp, r, t2, i)
      stats::dnorm(r_i, cm$mean, sqrt(cm$var + t2[i]), log = TRUE)
    }, numeric(1))
    elpd_i <- max(lp) + log(mean(exp(lp - max(lp))))
    pw <- lo$pointwise
    lpd_i <- pw[i, "elpd_loo"] + pw[i, "p_loo"]
    pw[i, "elpd_loo"] <- elpd_i
    pw[i, "p_loo"] <- lpd_i - elpd_i
    pw[i, "looic"] <- -2 * elpd_i
    if ("mcse_elpd_loo" %in% colnames(pw)) pw[i, "mcse_elpd_loo"] <- 0
    if ("influence_pareto_k" %in% colnames(pw)) pw[i, "influence_pareto_k"] <- 0
    lo$pointwise <- pw
    lo$diagnostics$pareto_k[i] <- 0
  }
  n <- nrow(lo$pointwise)
  cols <- c("elpd_loo", "p_loo", "looic")
  est <- t(vapply(cols, function(cl) {
    x <- lo$pointwise[, cl]
    c(Estimate = sum(x), SE = sqrt(n * stats::var(x)))
  }, numeric(2)))
  lo$estimates[cols, ] <- est
  for (cl in cols) lo[[cl]] <- est[cl, "Estimate"]
  for (cl in cols) lo[[paste0("se_", cl)]] <- est[cl, "SE"]
  attr(lo, "refitted") <- sp[idx]
  lo
}
