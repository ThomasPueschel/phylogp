# Prediction on the tree: posterior distributions at tips and internal nodes.
#
# Every branch-wise kernel, and the Ornstein-Uhlenbeck process, is a Gaussian
# Markov process along the tree: the value at a child node is
#   x_child = a_e * x_parent + eps_e,   eps_e ~ N(0, w_e),
# with a_e = 1 for branch-wise kernels and a_e = exp(-alpha * length) for OU.
# Observed tips add independent noise t2. One upward and one downward pass of
# Gaussian messages give the posterior of every node in O(N).

# Posterior moments of the latent process at every node.
# r: residuals at the tips (NA = not observed); a, w: per-edge coefficient and
# variance; rv: root variance (0 = root fixed at 0); t2: tip noise variances.
# Returns the posterior mean and variance at each node, and the "outside"
# moments (distribution of a node given the data outside its subtree).
tree_posterior <- function(ti, r, w, rv, t2, a = rep(1, ti$E)) {
  N <- ti$N; nn <- ti$n_nodes; E <- ti$E
  obs <- !is.na(r)
  exact <- obs & t2 <= 0                       # observed without noise
  P <- numeric(nn); H <- numeric(nn)           # information from the subtree
  P[seq_len(N)][obs & !exact] <- 1 / t2[obs & !exact]
  H[seq_len(N)][obs & !exact] <- r[obs & !exact] / t2[obs & !exact]
  msg_m <- numeric(E); msg_s <- numeric(E)     # message from child, as a function of the parent
  for (e in seq_len(E)) {
    c <- ti$des[e]; p <- ti$anc[e]
    if (c <= N && exact[c]) {
      m <- r[c]; s <- 0
    } else {
      s <- if (P[c] > 0) 1 / P[c] else Inf
      m <- if (P[c] > 0) H[c] / P[c] else 0
    }
    msg_m[e] <- m / a[e]
    msg_s[e] <- (s + w[e]) / a[e]^2
    if (is.finite(msg_s[e])) {
      P[p] <- P[p] + 1 / msg_s[e]
      H[p] <- H[p] + msg_m[e] / msg_s[e]
    }
  }
  om <- numeric(nn); ov <- numeric(nn)
  ov[ti$root] <- rv
  for (e in rev(seq_len(E))) {
    c <- ti$des[e]; p <- ti$anc[e]
    if (ov[p] <= 0) {
      mp <- om[p]; vp <- 0
    } else {
      own <- if (is.finite(msg_s[e])) c(1 / msg_s[e], msg_m[e] / msg_s[e]) else c(0, 0)
      Pp <- 1 / ov[p] + P[p] - own[1]
      Hp <- om[p] / ov[p] + H[p] - own[2]
      mp <- Hp / Pp; vp <- 1 / Pp
    }
    om[c] <- a[e] * mp
    ov[c] <- a[e]^2 * vp + w[e]
  }
  pm <- om; pv <- ov
  inner <- ov > 0
  prec <- 1 / ov[inner] + P[inner]
  pm[inner] <- (om[inner] / ov[inner] + H[inner]) / prec
  pv[inner] <- 1 / prec
  ex <- which(exact)
  pm[ex] <- r[ex]; pv[ex] <- 0
  list(mean = pm, var = pv, out_mean = om, out_var = ov)
}

# Moments of the latent process at tips `idx` given the observed residuals
# r (NA = unobserved): linear-time for branch-wise kernels, dense otherwise.
condition_tips <- function(ti, kp, r, t2, idx) {
  if (is.null(kp$K)) {
    post <- tree_posterior(ti, r, kp$w, kp$rv, t2)
    return(list(mean = post$mean[idx], var = post$var[idx]))
  }
  o <- which(!is.na(r))
  Kp <- edge_matrix(ti, kp$w, kp$rv, rep(0, ti$N)) + kp$K
  Koo <- Kp[o, o]
  diag(Koo) <- diag(Koo) + t2[o]
  U <- chol(Koo)
  A <- t(backsolve(U, backsolve(U, t(Kp[idx, o, drop = FALSE]), transpose = TRUE)))
  list(mean = as.numeric(A %*% r[o]),
       var = diag(Kp)[idx] - rowSums(A * Kp[idx, o, drop = FALSE]))
}

# Several observations per species: the within-species error of each species
# mean (added to the tip variances) and the log density of each species'
# deviations from its mean, log p(y_i | ybar_i), for one draw.
within_t2 <- function(fit, row) {
  if (is.null(fit$repeated)) return(0)
  unname(row["within_sigma"])^2 / fit$repeated$n
}
within_lp <- function(fit, row, beta) {
  rp <- fit$repeated
  if (is.null(rp)) return(rep(0, length(fit$tree$tip.label)))
  s2 <- unname(row["within_sigma"])^2
  d <- rp$y_dev - as.numeric(rp$X_dev %*% beta)
  ss <- as.numeric(rowsum(d^2, rp$species, reorder = TRUE))
  -0.5 * (rp$n - 1) * log(2 * pi * s2) - 0.5 * log(rp$n) - ss / (2 * s2)
}

# Parameter draws of a fit as a matrix (one row per draw)
fit_draws <- function(fit, ndraws = NULL, seed = 1) {
  if (fit$method == "optimize") {
    m <- fit$fit$mle()
    return(matrix(m, 1, dimnames = list(NULL, names(m))))
  }
  d <- posterior::as_draws_matrix(fit$fit$draws())
  if (!is.null(ndraws) && ndraws < nrow(d)) {
    set.seed(seed)
    d <- d[sort(sample(nrow(d), ndraws)), , drop = FALSE]
  }
  unclass(d)[, , drop = FALSE]
}

# Kernel parameters of one draw, as the list kernel_parts() expects
draw_pars <- function(L, row, ti) {
  take <- function(prefix, n) unname(row[paste0(prefix, "[", seq_len(n), "]")])
  out <- lapply(seq_along(L$ids), function(j) {
    cmp <- L$kernel$components[[j]]; id <- L$ids[j]
    p <- list(sigma = unname(row[paste0(id, "_sigma")]))
    for (nm in names(cmp$pars)) p[[nm]] <- unname(row[paste0(id, "_", nm)])
    if (identical(cmp$name, "relaxed")) {
      p$tau <- unname(row[paste0(id, "_tau")])
      p$z <- take(paste0(id, "_z"), ti$n_nodes)
    }
    if (identical(cmp$name, "ratetime")) p$lograte <- take(paste0(id, "_lograte"), cmp$n_bins)
    if (!is.null(cmp$draw_pars)) p <- cmp$draw_pars(p, row, id, ti)
    p
  })
  stats::setNames(out, L$ids)
}

# Edge coefficients for the OU process itself (rather than its tip kernel),
# used for ancestral states: only for a single OU component (+ white noise).
ou_transitions <- function(L, P, ti) {
  j <- which(vapply(L$kernel$components, `[[`, "", "name") == "ou")
  p <- P[[j]]; cmp <- L$kernel$components[[j]]
  a <- exp(-p$alpha * ti$len)
  list(a = a, w = p$sigma^2 / (2 * p$alpha) * (1 - a^2),
       rv = if (cmp$ou_root == "random") p$sigma^2 / (2 * p$alpha) else 0)
}

# Summarise per-draw Gaussian moments (rows = draws) into a data frame
summarise_moments <- function(M, V, probs, names, n_rep = 4) {
  S <- nrow(M)
  mean <- colMeans(M)
  sd <- sqrt(colMeans(V) + if (S > 1) apply(M, 2, function(x) mean((x - mean(x))^2)) else 0)
  out <- data.frame(name = names, mean = mean, sd = sd)
  if (S == 1) {
    q <- vapply(probs, function(pr) stats::qnorm(pr, M[1, ], sqrt(V[1, ])), numeric(ncol(M)))
  } else {
    sims <- do.call(rbind, lapply(seq_len(n_rep), function(k)
      M + sqrt(V) * matrix(stats::rnorm(length(M)), nrow(M))))
    q <- t(apply(sims, 2, stats::quantile, probs = probs, names = FALSE))
  }
  q <- matrix(q, ncol = length(probs))
  for (k in seq_along(probs)) out[[sprintf("q%g", 100 * probs[k])]] <- q[, k]
  out
}

#' Predict trait values for species
#'
#' Posterior predictive distributions of the trait for species in the tree:
#' species that were fitted (smoothed values) or species in the tree without
#' trait data (phylogenetic imputation, i.e. kriging on the tree). Each
#' prediction conditions on all observed species, through the kernel and the
#' regression part of the model, and integrates over the posterior draws.
#'
#' @param object A `phylogp_fit`.
#' @param newdata Optional data frame with the predictor values of the species
#'   to predict, named as in the fit (row names or the `species` column). Its
#'   species must be in `tree`. If `NULL`, the fitted species are predicted.
#' @param tree Tree containing the species to predict; by default the tree
#'   given to [phylo_gp()].
#' @param type `"response"` (a new observation of the trait, including
#'   residual variation from `k_white()` and, for fits with several rows per
#'   species, the within-species variation) or `"latent"` (the evolutionary
#'   process at the tip).
#' @param ndraws Number of posterior draws used.
#' @param probs Quantiles to report.
#' @param seed Random seed for draw selection and quantiles.
#' @param ... Not used.
#' @return A data frame with one row per species: whether it was observed,
#'   posterior mean, standard deviation and quantiles.
#' @examples
#' \dontrun{
#' fit <- phylo_gp(y ~ x, d, tree, k_bm() + k_white())
#' predict(fit, newdata = data.frame(x = 1.2, row.names = "missing_species"))
#' }
#' @references
#' Garland T. Jr., Ives A.R. 2000. Using the past to predict the present: Confidence intervals for regression equations in phylogenetic comparative methods. The American Naturalist. 155:346–364.
#'
#' Rasmussen C.E., Williams C.K.I. 2006. Gaussian processes for machine learning. Cambridge, MA: MIT Press.
#' @export
predict.phylogp_fit <- function(object, newdata = NULL, tree = NULL,
                                type = c("response", "latent"), ndraws = 500,
                                probs = c(0.05, 0.95), seed = 1, ...) {
  type <- match.arg(type)
  full <- tree %||% object$tree_full %||% object$tree
  obs_sp <- object$tree$tip.label
  y_obs <- stats::setNames(object$standata$y, obs_sp)
  se_obs <- stats::setNames(object$standata$se, obs_sp)
  X_obs <- object$standata$X
  if (is.null(newdata)) {
    targets <- obs_sp
    X_t <- X_obs
  } else {
    newdata <- as.data.frame(newdata)
    targets <- if (is.null(object$species)) rownames(newdata) else as.character(newdata[[object$species]])
    miss <- setdiff(targets, full$tip.label)
    if (length(miss)) stop("Species not in the tree: ", paste(utils::head(miss, 3), collapse = ", "), call. = FALSE)
    tt <- stats::delete.response(object$terms)
    mf <- stats::model.frame(tt, newdata, xlev = object$xlevels, na.action = stats::na.pass)
    X_t <- stats::model.matrix(tt, mf, xlev = object$xlevels)
    if (anyNA(X_t)) stop("Missing predictor values in `newdata`.", call. = FALSE)
  }
  new_sp <- setdiff(targets, obs_sp)
  L <- object$layout
  if (length(new_sp) && L$latent)
    stop("Predicting species that were not fitted is not yet available for kernels ",
         "with latent values (k_relaxed(), k_rate_time(), latent slopes, k_white(df = )).",
         call. = FALSE)
  kernel <- object$kernel
  if (length(new_sp) && !is.null(newdata)) kernel <- extend_slopes(kernel, newdata, targets)
  ti <- if (length(new_sp)) {
    keep <- union(obs_sp, new_sp)
    tree_info(ape::drop.tip(full, setdiff(full$tip.label, keep)))
  } else tree_info(object$tree)
  if (length(new_sp) && abs(ti$Tmax - tree_info(object$tree)$Tmax) > 1e-8 * ti$Tmax)
    warning("Adding the new species changes the tree height; kernels that use it ",
            "(delta, early burst, OU) are evaluated on the extended tree.", call. = FALSE)
  L <- kernel_layout(kernel, ti)
  o <- match(obs_sp, ti$labels)
  t_idx <- match(targets, ti$labels)
  D <- fit_draws(object, ndraws, seed)
  ncoef <- length(object$coef_names)
  M <- V <- matrix(0, nrow(D), length(targets))
  for (s in seq_len(nrow(D))) {
    row <- D[s, ]
    beta <- unname(row[paste0("coefs[", seq_len(ncoef), "]")])
    kp <- kernel_parts(L$kernel, ti, draw_pars(L, row, ti))
    r <- rep(NA_real_, ti$N)
    r[o] <- y_obs - as.numeric(X_obs %*% beta) - kp$shift[o]
    t2 <- kp$t2
    t2[o] <- t2[o] + se_obs^2 + within_t2(object, row)
    cm <- condition_tips(ti, kp, r, t2, t_idx)
    m <- cm$mean + kp$shift[t_idx]; v <- cm$var
    if (type == "response")
      v <- v + kp$t2[t_idx] + if (is.null(object$repeated)) 0 else unname(row["within_sigma"])^2
    M[s, ] <- m + as.numeric(X_t %*% beta)
    V[s, ] <- pmax(v, 0)
  }
  set.seed(seed)
  out <- summarise_moments(M, V, probs, targets)
  names(out)[1] <- "species"
  out <- cbind(out[1], observed = targets %in% obs_sp, out[-1])
  rownames(out) <- NULL
  out
}

#' Ancestral states
#'
#' Posterior distributions of the trait at the internal nodes of the tree
#' (kriging at internal points of the Gaussian process), computed in linear
#' time and averaged over posterior draws. Available for intercept-only models
#' whose kernel is a process along the whole tree: sums of branch-wise
#' components (plus residual or measurement error at the tips), or a single
#' Ornstein-Uhlenbeck component. For Pagel's kappa and delta, reconstructions
#' refer to the transformed tree; for Pagel's lambda, to the tree with
#' internal branches scaled by lambda.
#'
#' For OU, reconstructions use the OU process itself, not its branch-wise tip
#' kernel: on an ultrametric tree OU and an accelerating early burst have the
#' same tip distribution but different ancestral states.
#'
#' @param object A `phylogp_fit` of an intercept-only model (`y ~ 1`).
#' @param ndraws Number of posterior draws used.
#' @param probs Quantiles to report.
#' @param tips Also return the tips?
#' @param seed Random seed.
#' @return A data frame with one row per node (numbered as in `object$tree`):
#'   posterior mean, standard deviation and quantiles.
#' @references
#' Schluter D., Price T., Mooers A.Ø., Ludwig D. 1997. Likelihood of ancestor states in adaptive radiation. Evolution. 51:1699–1711.
#'
#' Garland T. Jr., Ives A.R. 2000. Using the past to predict the present: Confidence intervals for regression equations in phylogenetic comparative methods. The American Naturalist. 155:346–364.
#' @export
ancestral <- function(object, ndraws = 500, probs = c(0.05, 0.95), tips = FALSE, seed = 1) {
  if (!inherits(object, "phylogp_fit")) stop("`object` must be a phylogp_fit.", call. = FALSE)
  if (!identical(object$coef_names, "(Intercept)"))
    stop("Ancestral states need an intercept-only model (y ~ 1).", call. = FALSE)
  ti <- tree_info(object$tree)
  L <- kernel_layout(object$kernel, ti)
  nm <- vapply(L$kernel$components, `[[`, "", "name")
  is_ou <- nm == "ou"
  if (any(is_ou) && !all(nm[!is_ou] == "white"))
    stop("Ancestral states with OU are available for k_ou() alone or k_ou() + k_white().", call. = FALSE)
  if (any(L$forms == "slope"))
    stop("Ancestral states are not available for models with varying slopes.", call. = FALSE)
  if (!any(is_ou) && L$dense)
    stop("Ancestral states need a kernel that is a process along the tree ",
         "(no dense components such as k_dist()).", call. = FALSE)
  D <- fit_draws(object, ndraws, seed)
  y <- object$standata$y; se <- object$standata$se
  nodes <- if (tips) seq_len(ti$n_nodes) else (ti$N + 1):ti$n_nodes
  M <- V <- matrix(0, nrow(D), length(nodes))
  for (s in seq_len(nrow(D))) {
    row <- D[s, ]
    mu <- unname(row["coefs[1]"])
    P <- draw_pars(L, row, ti)
    kp <- kernel_parts(L$kernel, ti, P)
    t2 <- kp$t2 + se^2 + within_t2(object, row)
    post <- if (any(is_ou)) {
      tr <- ou_transitions(L, P, ti)
      tree_posterior(ti, y - mu, tr$w, tr$rv, t2, a = tr$a)
    } else {
      tree_posterior(ti, y - mu, kp$w, kp$rv, t2)
    }
    M[s, ] <- mu + post$mean[nodes]
    V[s, ] <- pmax(post$var[nodes], 0)
  }
  set.seed(seed)
  out <- summarise_moments(M, V, probs, nodes)
  names(out)[1] <- "node"
  out
}
