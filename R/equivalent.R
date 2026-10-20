# When do two kernels imply the same distribution for the tips?

#' Do two kernels imply the same distribution for the tips?
#'
#' Different evolutionary processes can imply exactly, or nearly, the same
#' distribution for the species at the tips of a tree, and then no amount of
#' extant data can tell them apart. On an ultrametric tree, for example, an
#' Ornstein-Uhlenbeck process with a fixed root has the same kernel as an
#' early-burst model whose accelerating rate is twice the OU `alpha`.
#'
#' For each of several parameter settings of `k1`, `kernel_equivalent()`
#' searches over all parameters of `k2` (including the scales `sigma`) for the
#' tip distribution closest to that of `k1`, measured by the Kullback-Leibler
#' divergence
#' \deqn{KL = E_{k1}[\log p_{k1}(y) - \log p_{k2}(y)].}
#' KL is the expected log-likelihood ratio in favour of `k1` for one data set
#' of these tips simulated from `k1`, when `k2`'s parameters are chosen as well
#' as possible. Values near 0 mean `k2` can reproduce `k1`; values below about
#' 1 mean the two will rarely be distinguishable with data of this size.
#'
#' With `mean = "free"` (default), the trait mean is an unknown parameter, so
#' only contrasts between species are compared, as in REML. Two kernels that
#' differ only by a constant added to every covariance (such as a random
#' versus a fixed root) are then equivalent. With `mean = "zero"` the full
#' covariance matrices are compared.
#'
#' The search is numerical, with several starting points. A small KL is
#' therefore evidence of equivalence, but a large one means only that no
#' equivalent parameter values were found.
#'
#' @param k1 Kernel whose tip distribution is to be reproduced.
#' @param k2 Candidate kernel.
#' @param tree A `phylo` tree (or the output of [tree_info()]).
#' @param pars Optional parameter values of `k1` (as in [kernel_matrix()]); by
#'   default a grid of `n_grid` settings spanning the bounds of its parameters.
#' @param n_grid Number of parameter settings of `k1` tried.
#' @param mean `"free"` (unknown mean, contrasts only) or `"zero"`.
#' @param n_starts Number of starting points of the search over `k2`.
#' @param seed Random seed for the starting points.
#' @return A data frame (class `phylogp_equivalence`) with one row per setting
#'   of `k1`: its parameters, the smallest KL found, the best parameters of
#'   `k2`, and a verdict (`"identical"`, `"nearly identical"` or
#'   `"distinguishable"`).
#' @examples
#' tr <- ape::rcoal(30)
#' # OU with a fixed root vs accelerating early burst: identical on this tree
#' kernel_equivalent(k_ou(root = "fixed"), k_eb(upper = 20), tr, n_grid = 3)
#' @export
kernel_equivalent <- function(k1, k2, tree, pars = NULL, n_grid = 8,
                              mean = c("free", "zero"), n_starts = 5, seed = 1) {
  mean <- match.arg(mean)
  ti <- tree_info(tree)
  L1 <- kernel_layout(k1, ti)
  L2 <- kernel_layout(k2, ti)
  if (L1$latent || L2$latent)
    stop("Kernels with latent values (relaxed rates, rate through time, latent slopes, ",
         "Student-t residuals) have no fixed tip covariance to compare.", call. = FALSE)
  Q <- if (mean == "free") qr.Q(qr(matrix(1, ti$N, 1)), complete = TRUE)[, -1, drop = FALSE] else NULL
  project <- function(K) if (is.null(Q)) K else crossprod(Q, K %*% Q)
  cov_of <- function(L, P) {
    kp <- kernel_parts(L$kernel, ti, P)
    K <- edge_matrix(ti, kp$w, kp$rv, kp$t2)
    if (!is.null(kp$K)) K <- K + kp$K
    project(K)
  }

  settings <- if (!is.null(pars)) list(pars) else k1_grid(L1, ti, n_grid)
  space <- k2_space(L2, ti)
  set.seed(seed)
  rows <- lapply(settings, function(P1) {
    A <- cov_of(L1, P1)
    UA <- tryCatch(chol(A), error = function(e) NULL)
    if (is.null(UA)) return(NULL)                      # k1 not positive definite here
    logdetA <- 2 * sum(log(diag(UA)))
    m <- nrow(A)
    kl <- function(phi) {
      B <- cov_of(L2, space$to_pars(phi))
      UB <- tryCatch(chol(B), error = function(e) NULL)
      if (is.null(UB)) return(1e10)
      X <- backsolve(UB, t(UA), transpose = TRUE)
      0.5 * (sum(X^2) - m + 2 * sum(log(diag(UB))) - logdetA)
    }
    # starting points: k2 defaults, rescaled to the overall variance of k1
    base <- space$start
    s <- mean(diag(A)) / mean(diag(cov_of(L2, space$to_pars(base))))
    base[space$is_scale] <- base[space$is_scale] + 0.5 * log(s)
    starts <- c(list(base), lapply(seq_len(n_starts - 1), function(k)
      base + stats::rnorm(length(base), 0, 1.5)))
    fits <- lapply(starts, function(st) tryCatch(
      stats::optim(st, kl, method = "BFGS", control = list(reltol = 1e-14, maxit = 2000)),
      error = function(e) NULL))
    fits <- Filter(Negate(is.null), fits)
    best <- fits[[which.min(vapply(fits, `[[`, 0, "value"))]]
    q <- stats::plogis(best$par[!space$is_scale])
    data.frame(as.list(flatten_pars(L1, P1, "k1")),
               KL = max(best$value, 0),
               as.list(flatten_pars(L2, space$to_pars(best$par), "k2")),
               k2_at_bound = any(q < 1e-4 | q > 1 - 1e-4),
               check.names = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  out$verdict <- ifelse(out$KL < 1e-6, "identical",
                        ifelse(out$KL < 1, "nearly identical", "distinguishable"))
  structure(out, class = c("phylogp_equivalence", "data.frame"),
            labels = c(kernel_label(L1$kernel), kernel_label(L2$kernel)),
            mean = mean, n_tips = ti$N, ultrametric = ti$ultrametric)
}

kernel_label <- function(kernel)
  paste(vapply(kernel$components, `[[`, "", "label"), collapse = " + ")

# Grid of parameter settings for k1, varying all components together.
# Scale-type parameters (bounds spanning orders of magnitude) are spread over
# the top three orders of magnitude of their range, where kernels differ.
k1_grid <- function(L, ti, n) {
  grids <- lapply(L$kernel$components, function(cmp) {
    for (p in names(cmp$pars)) {
      b <- resolve_bounds(cmp$pars[[p]], ti, cmp)
      if (b[1] > 0 && b[2] / b[1] > 1e3) cmp$pars[[p]] <- bounded(b[2] * 1e-3, b[2])
    }
    par_grid(cmp, ti, n)
  })
  n_set <- max(lengths(grids))
  lapply(seq_len(n_set), function(g)
    stats::setNames(lapply(grids, function(gr) gr[[(g - 1) %% length(gr) + 1]]), L$ids))
}

# Unconstrained parameterisation of k2: log sigma for every component,
# logit-scaled bounded parameters
k2_space <- function(L, ti) {
  entries <- list()
  for (j in seq_along(L$ids)) {
    cmp <- L$kernel$components[[j]]
    entries[[length(entries) + 1]] <- list(j = j, name = "sigma", scale = TRUE)
    for (p in names(cmp$pars))
      entries[[length(entries) + 1]] <- list(j = j, name = p, scale = FALSE,
                                             b = resolve_bounds(cmp$pars[[p]], ti, cmp))
  }
  defaults <- lapply(L$kernel$components, function(cmp) cmp$defaults(ti))
  start <- vapply(entries, function(en) {
    v <- defaults[[en$j]][[en$name]]
    if (en$scale) return(log(v))
    if (is.null(v) || is.na(v)) v <- mean(en$b)
    q <- (v - en$b[1]) / diff(en$b)
    stats::qlogis(min(max(q, 0.01), 0.99))
  }, 0)
  to_pars <- function(phi) {
    P <- stats::setNames(lapply(defaults, function(d) d), L$ids)
    for (k in seq_along(entries)) {
      en <- entries[[k]]
      P[[en$j]][[en$name]] <- if (en$scale) exp(phi[k]) else
        en$b[1] + diff(en$b) * stats::plogis(phi[k])
    }
    P
  }
  list(start = start, to_pars = to_pars,
       is_scale = vapply(entries, `[[`, NA, "scale"))
}

# Named vector of the parameter values of a kernel, e.g. k1:ou1.alpha
flatten_pars <- function(L, P, prefix) {
  out <- c()
  for (j in seq_along(L$ids)) {
    cmp <- L$kernel$components[[j]]
    for (nm in c("sigma", names(cmp$pars))) {
      v <- P[[L$ids[j]]][[nm]]
      if (is.numeric(v) && length(v) == 1)
        out[paste0(prefix, ":", L$ids[j], ".", nm)] <- v
    }
  }
  out
}

#' @export
print.phylogp_equivalence <- function(x, digits = 3, ...) {
  lab <- attr(x, "labels")
  cat(sprintf("Can %s reproduce the tip distribution of %s?\n", lab[2], lab[1]))
  cat(sprintf("%d tips, %s tree, %s mean\n\n", attr(x, "n_tips"),
              if (attr(x, "ultrametric")) "ultrametric" else "non-ultrametric",
              if (attr(x, "mean") == "free") "unknown (contrasts compared)" else "zero"))
  df <- as.data.frame(unclass(x))
  attr(df, "labels") <- attr(df, "mean") <- attr(df, "n_tips") <- attr(df, "ultrametric") <- NULL
  num <- vapply(df, is.numeric, NA)
  df[num] <- lapply(df[num], signif, digits)
  print(df, row.names = FALSE)
  cat("\nKL: expected log-likelihood ratio in favour of the first kernel for one\n",
      "data set simulated from it, with the second kernel's parameters chosen as\n",
      "well as possible.\n", sep = "")
  v <- unique(x$verdict)
  if (any(x$k2_at_bound & x$verdict != "identical"))
    cat("In some settings the best parameters of the second kernel are at their bounds;\n",
        "widen the bounds (e.g. k_eb(upper = )) before concluding that the kernels differ.\n", sep = "")
  msg <- if (identical(v, "identical")) {
    "Identical for every parameter setting tried: the tips cannot distinguish the two kernels."
  } else if (all(v %in% c("identical", "nearly identical"))) {
    "Nearly identical for every setting tried: expect weak evidence either way."
  } else if (all(v == "distinguishable")) {
    "No equivalent parameter values were found for any setting tried."
  } else {
    "Equivalence depends on the parameter values of the first kernel."
  }
  cat(msg, "\n")
  invisible(x)
}
