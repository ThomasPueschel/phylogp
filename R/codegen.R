# Stan code generation.
#
# The generated program has one structure for every kernel:
#   * edge (branch-wise) and white components add to the edge variances w,
#     the root variance rv and the tip variances t2;
#   * if every component is of this kind the likelihood is the O(N) pruning
#     recursion and leave-one-out densities come from an O(N) two-pass
#     message passing; otherwise the dense kernel K = B diag(w) B' + rv + t2
#     plus the dense components is used (O(N^3)).

sub_id <- function(x, id) gsub("ID_", paste0(id, "_"), x, fixed = TRUE)

# Forms of the components on a given tree (OU is branch-wise only on
# ultrametric trees). Without a tree, assume an ultrametric one.
component_forms <- function(kernel, ti = NULL) {
  probe <- ti %||% list(ultrametric = TRUE)
  vapply(kernel$components, function(cmp) cmp$form(probe), "")
}

kernel_layout <- function(kernel, ti = NULL) {
  kernel <- as_kernel(kernel)
  forms <- component_forms(kernel, ti)
  ids <- component_ids(kernel)
  kinds <- vapply(kernel$components, function(cmp) cmp$tip_loo$kind %||% "", "")
  quad <- which(kinds %in% c("gh", "lam"))
  if (length(quad) > 1)
    stop("At most one component whose leave-one-out needs quadrature ",
         "(k_relaxed(), k_white(df = )) is supported.", call. = FALSE)
  dense <- any(forms == "dense")
  needs <- unique(unlist(lapply(kernel$components[forms == "dense"], `[[`, "needs")))
  list(kernel = kernel, forms = forms, ids = ids, dense = dense, needs = needs,
       quad = quad, quad_kind = if (length(quad)) kinds[quad] else "",
       closed = which(kinds == "closed"),
       latent = any(vapply(kernel$components, `[[`, NA, "latent")))
}

#' Stan code for a kernel or a fitted model
#'
#' Writes the complete Stan program that [phylo_gp()] compiles for a kernel.
#' The program can be inspected, edited and run directly with `cmdstanr`
#' together with the data from [kernel_standata()]. Given a fitted model, it
#' returns the program that was fitted.
#'
#' The name differs from the `make_stancode` and `stancode` functions of
#' brms, so that both packages can be loaded together.
#'
#' @param kernel A kernel, e.g. `k_bm() + k_white()`, or a `phylogp_fit`.
#' @param tree Optional tree. It decides whether components such as
#'   `k_ou()` can use the branch-wise form (ultrametric trees) or need dense
#'   algebra. Without a tree an ultrametric tree is assumed. Ignored for a
#'   fitted model.
#' @param repeated Write the program for several observations per species?
#'   The kernel then describes the species means, and a within-species
#'   residual standard deviation `within_sigma` is added; see [phylo_gp()].
#' @return A character string with class `phylogp_stancode`.
#' @examples
#' cat(kernel_stancode(k_bm() + k_white()))
#' @export
kernel_stancode <- function(kernel, tree = NULL, repeated = FALSE) {
  if (inherits(kernel, "phylogp_fit")) return(kernel$code)
  ti <- if (is.null(tree)) NULL else tree_info(tree)
  L <- kernel_layout(kernel, ti)
  comps <- L$kernel$components
  ind <- function(x, n = 1) if (length(x)) paste0(strrep("  ", n), x) else character()

  # ---- data ----
  data <- c(
    "int<lower=1> N;", "int<lower=1> P;", "matrix[N, P] X;", "vector[N] y;",
    "vector<lower=0>[N] se;                // known measurement error SDs",
    "// tree, edges in postorder (children before parents)",
    "int<lower=1> E;", "int<lower=N + 1> n_nodes;", "int<lower=N + 1> root;",
    "array[E] int anc;", "array[E] int des;", "array[E] int first_child;",
    "array[N] int tip_edge;",
    "vector<lower=0>[E] len;", "vector<lower=0>[E] hs;", "vector<lower=0>[E] he;",
    "vector<lower=0, upper=1>[E] is_tip;", "real<lower=0> Tmax;",
    if (L$dense) "matrix[N, E] B;                       // tip-by-edge incidence",
    if ("tC" %in% L$needs) "matrix[N, N] tC;                      // shared path lengths t_ij",
    if ("Dtree" %in% L$needs) "matrix[N, N] Dtree;                   // patristic distances d_ij",
    "int<lower=0, upper=1> use_prior;      // 0: posterior mode = maximum likelihood",
    "real<lower=0> prior_beta_scale;",
    if (L$quad_kind == "gh") c("int<lower=1> n_gh;", "vector[n_gh] gh_x;", "vector[n_gh] gh_w;"),
    if (L$quad_kind == "lam") c("int<lower=1> n_lq;", "vector<lower=0>[n_lq] lq_lam;", "vector<lower=0>[n_lq] lq_w;"),
    if (repeated) c(
      "// several observations per species (y and X hold species means)",
      "int<lower=1> M;                       // observations",
      "array[M] int<lower=1, upper=N> obs_sp; // species of each observation",
      "vector[M] y_obs;                      // response of each observation",
      "matrix[M, P] X_obs;                   // design matrix of the observations",
      "vector<lower=0>[M] se_obs;            // known measurement error of each observation",
      "vector<lower=0, upper=1>[N] within_include; // 0 leaves a species out",
      "real<lower=0> within_sigma_scale;"))
  for (j in seq_along(comps)) {
    cmp <- comps[[j]]; id <- L$ids[j]
    data <- c(data, sprintf("// %s: %s", id, cmp$label),
              sprintf("real<lower=0> %s_sigma_scale;", id))
    for (p in names(cmp$pars))
      data <- c(data, sprintf("real %1$s_%2$s_lb;", id, p),
                sprintf("real<lower=%1$s_%2$s_lb> %1$s_%2$s_ub;", id, p))
    data <- c(data, sub_id(cmp$stan_data, id),
              if (repeated && L$forms[j] == "slope")
                sprintf("vector[M] %s_xs_obs;                // x of each observation, centred", id))
  }

  # ---- parameters ----
  pars <- c("vector[P] coefs;              // regression coefficients",
            if (repeated) "real<lower=0> within_sigma;    // residual sd within species")
  tpars <- character()
  for (j in seq_along(comps)) {
    cmp <- comps[[j]]; id <- L$ids[j]
    pars <- c(pars, sprintf("real<lower=0> %s_sigma;", id))
    for (p in names(cmp$pars))
      pars <- c(pars, sprintf("real<lower=%1$s_%2$s_lb, upper=%1$s_%2$s_ub> %1$s_%2$s;", id, p))
    pars <- c(pars, sub_id(cmp$stan_pars, id))
    tpars <- c(tpars, sub_id(cmp$stan_tpars, id))
  }

  # ---- kernel assembly (used in model and generated quantities) ----
  asm <- c("vector[E] w = rep_vector(0, E);", "real rv = 0;", "vector[N] t2 = square(se);",
           if (repeated) c(
             "vector[M] r_obs = y_obs - X_obs * coefs;",
             sprintf("r_obs -= %1$s_xs_obs .* %1$s_v[obs_sp];  // varying slope of each row",
                     L$ids[L$forms == "slope"]),
             "// species means of the residuals, weighted by 1 / (within_sigma^2 + se_obs^2),",
             "// their error variances and the log density of the rows given the means",
             "matrix[N, 3] wp = within_parts(r_obs, se_obs, obs_sp, N, within_sigma,",
             "                               rep_vector(0, M), rep_vector(0, N));",
             "t2 += wp[, 2];"))
  for (j in seq_along(comps)) {
    cmp <- comps[[j]]; id <- L$ids[j]
    if (L$forms[j] == "white") {
      asm <- c(asm, if (is.null(cmp$white_expr)) sprintf("t2 += square(%s_sigma);", id) else
        sprintf("t2 += square(%s_sigma) * (%s);", id, sub_id(cmp$white_expr, id)))
    } else if (L$forms[j] == "edge") {
      code <- sprintf("w += square(%s_sigma) * (%s);", id, sub_id(cmp$edge, id))
      if (!is.null(cmp$root))
        code <- c(code, sprintf("rv += square(%s_sigma) * (%s);", id, sub_id(cmp$root, id)))
      asm <- c(asm, sprintf("// %s", cmp$label), with_aliases(id, cmp, code))
    }
  }
  if (L$dense) {
    asm <- c(asm, "matrix[N, N] K = tree_kernel(B, w, rv, t2);")
    for (j in which(L$forms == "dense")) {
      cmp <- comps[[j]]; id <- L$ids[j]
      code <- sprintf("K += square(%s_sigma) * (%s);", id, sub_id(cmp$dense, id))
      asm <- c(asm, sprintf("// %s", cmp$label),
               with_aliases(id, cmp, code, extra = sub_id(cmp$dense_alias, id)))
    }
  }
  asm <- c(asm, if (repeated) "vector[N] r = wp[, 1];" else c("vector[N] r = y - X * coefs;",
           unlist(lapply(seq_along(comps), function(j) sub_id(comps[[j]]$stan_shift, L$ids[j])))))

  # ---- model ----
  lik <- if (L$dense)
    "target += multi_normal_cholesky_lpdf(r | rep_vector(0, N), cholesky_decompose(K));" else
    "target += tree_pruning_lpdf(r | w, rv, t2, anc, des, first_child, n_nodes, root);"
  latent_prior <- unlist(lapply(seq_along(comps), function(j)
    sub_id(comps[[j]]$stan_latent_prior, L$ids[j])))
  priors <- c("coefs ~ normal(0, prior_beta_scale);",
              unlist(lapply(seq_along(comps), function(j) {
                id <- L$ids[j]
                c(sprintf("%1$s_sigma ~ normal(0, %1$s_sigma_scale);", id),
                  sub_id(comps[[j]]$stan_prior, id))
              })))
  within <- "wp[, 3]"
  if (repeated) {
    lik <- c(lik, sprintf("target += dot_product(within_include, %s);", within))
    priors <- c(priors, "within_sigma ~ normal(0, within_sigma_scale);")
  }
  model <- c("{", ind(c(asm, lik)), "}", latent_prior,
             "if (use_prior == 1) {", ind(priors), "}")

  # ---- generated quantities: leave-one-out densities p(y_i | y_-i, theta) ----
  # mo: mean and variance of r_i given the other species. A left-out species'
  # own latent values (its slope or rate on the tip edge) are informed only by
  # y_i, so they are integrated out against their conditional prior: in closed
  # form (latent slopes) or by quadrature (relaxed rates, heavy-tailed residuals).
  loo_head <- if (L$dense) "matrix[N, 2] mo = gaussian_loo_moments(r, K);" else
    "matrix[N, 2] mo = tree_loo_moments(r, w, rv, t2, anc, des, n_nodes, root);"
  body <- c("int e = tip_edge[i];", "real ri = r[i];",
            if (L$dense) "real vi = mo[i, 2];" else "real vi = mo[i, 2] + t2[i];")
  # several rows per species with a varying slope: the left-out species' rows
  # are scored jointly, with its slope at the tip integrated out (wl)
  wl <- character()
  if (repeated && length(L$closed)) {
    id <- L$ids[L$closed]
    wl <- c(sprintf("matrix[N, 3] wl = within_parts(r_obs + %1$s_xs_obs .* (%1$s_v[obs_sp] - %1$s_v[anc[tip_edge[obs_sp]]]),", id),
            sprintf("                               se_obs, obs_sp, N, within_sigma, %1$s_xs_obs,", id),
            sprintf("                               square(%1$s_sigma) * %1$s_wb[tip_edge]);", id))
    body <- c(body, "ri = wl[i, 1];", "vi += wl[i, 2] - wp[i, 2];")
    within <- "wl[, 3]"
  }
  for (j in if (length(wl)) integer() else L$closed) {
    tl <- comps[[j]]$tip_loo; id <- L$ids[j]
    body <- c(body, sprintf("ri += %s;", sub_id(tl$r, id)), sprintf("vi += %s;", sub_id(tl$v, id)))
  }
  if (length(L$quad)) {
    j <- L$quad; tl <- comps[[j]]$tip_loo; id <- L$ids[j]
    ng <- if (L$quad_kind == "gh") "n_gh" else "n_lq"
    wt <- if (L$quad_kind == "gh") "log(gh_w[g]) - 0.5 * log(pi())" else "log(lq_w[g])"
    body <- c(body, sub_id(tl$setup, id),
              sprintf("vector[%s] lp;", ng),
              sprintf("vi -= %s;", sub_id(tl$cur, id)),
              sprintf("for (g in 1:%s)", ng),
              sprintf("  lp[g] = normal_lpdf(ri | mo[i, 1], sqrt(vi + %s)) + %s;",
                      sub_id(tl$node, id), wt),
              "log_lik[i] = log_sum_exp(lp);")
  } else {
    body <- c(body, "log_lik[i] = normal_lpdf(ri | mo[i, 1], sqrt(vi));")
  }
  loo <- c(loo_head, wl, "for (i in 1:N) {", ind(body), "}",
           if (repeated) sprintf("log_lik += %s;  // leave out all rows of a species", within))
  gq <- c("vector[N] log_lik;", "{", ind(c(asm, loo)), "}")

  fun_file <- system.file("stan", "phylogp_functions.stan", package = "phylogp")
  funs <- readLines(fun_file, warn = FALSE)
  header <- c(
    "// Generated by the phylogp R package (phylogenetic Gaussian processes)",
    sprintf("// Kernel: %s", paste(vapply(comps, `[[`, "", "label"), collapse = " + ")),
    sprintf("// Likelihood: %s", if (L$dense) "dense (O(N^3))" else "pruning (O(N))"))
  code <- c(header,
            "functions {", funs, if (repeated) ind(within_lp_stan()), "}",
            "data {", ind(data), "}",
            "parameters {", ind(pars), "}",
            if (length(tpars)) c("transformed parameters {", ind(tpars), "}"),
            "model {", ind(model), "}",
            "generated quantities {", ind(gq), "}")
  structure(paste(code, collapse = "\n"), class = c("phylogp_stancode", "character"))
}

#' @export
print.phylogp_stancode <- function(x, ...) {
  cat(x, "\n")
  invisible(x)
}

# Stan function for several observations per species. The residuals of
# species i are r_ij = g_i + xs_ij u_i + e_ij, with g_i the species value,
# e_ij ~ N(0, s^2 + se_ij^2) and, optionally, a slope u_i ~ N(0, s2u_i)
# (s2u = 0: none). The likelihood of the rows given g_i factorises into
# N(rbar_i | g_i, 1 / P_i), with rbar_i the generalised least-squares mean
# (Sherman-Morrison for the rank-one slope term), and the rows given rbar_i.
# Returns, per species, rbar_i, 1 / P_i and log p(r_i | rbar_i). An R copy is
# within_parts() in predict.R.
within_lp_stan <- function() c(
  "matrix within_parts(vector r, vector se_obs, array[] int obs_sp, int N, real s,",
  "                    vector xs, vector s2u) {",
  "  int M = rows(r);",
  "  vector[M] prec = inv(square(s) + square(se_obs));",
  "  vector[N] A = rep_vector(0, N);    // sum prec",
  "  vector[N] B = rep_vector(0, N);    // sum prec x",
  "  vector[N] C = rep_vector(0, N);    // sum prec x^2",
  "  vector[N] Hr = rep_vector(0, N);   // sum prec r",
  "  vector[N] Hx = rep_vector(0, N);   // sum prec x r",
  "  vector[N] lp = rep_vector(0, N);",
  "  vector[N] Q = rep_vector(0, N);",
  "  vector[N] n = rep_vector(0, N);",
  "  matrix[N, 3] out;",
  "  for (m in 1:M) {",
  "    int i = obs_sp[m];",
  "    A[i] += prec[m]; B[i] += prec[m] * xs[m]; C[i] += prec[m] * square(xs[m]);",
  "    Hr[i] += prec[m] * r[m]; Hx[i] += prec[m] * xs[m] * r[m];",
  "    lp[i] += 0.5 * log(prec[m]); n[i] += 1;",
  "  }",
  "  vector[N] c = 1 + s2u .* C;",
  "  vector[N] P = A - s2u .* square(B) ./ c;",
  "  out[, 1] = (Hr - s2u .* B .* Hx ./ c) ./ P;",
  "  out[, 2] = inv(P);",
  "  for (m in 1:M) Q[obs_sp[m]] += prec[m] * square(r[m] - out[obs_sp[m], 1]);",
  "  Q -= s2u .* square(Hx - out[, 1] .* B) ./ c;",
  "  out[, 3] = lp - 0.5 * log(c) - 0.5 * Q - 0.5 * (n - 1) * log(2 * pi()) - 0.5 * log(P);",
  "  return out;",
  "}")

#' Data for the Stan program of a kernel
#'
#' The data list for the program written by [kernel_stancode()].
#'
#' @param kernel A kernel.
#' @param tree A `phylo` tree, or the output of [tree_info()].
#' @param y Response vector in the tip order of `tree`.
#' @param X Design matrix in the same order (default: intercept only).
#' @param se Known measurement-error standard deviations (default 0).
#' @param use_prior `0` for maximum likelihood (no priors), `1` for Bayesian.
#' @param prior_sigma_scale Optional scale of the half-normal priors on the
#'   component scales; by default `2 sd(y) / sqrt(base variance)`.
#' @param prior_beta_scale Optional scale of the normal priors on the
#'   regression coefficients.
#' @param repeated For several observations per species: `NULL`, or a list
#'   with `n` (observations per species, in tip order), `species` (tip index
#'   of each observation), `y` and `X` (response and design matrix of the
#'   observations), optionally `se` (known measurement error of each
#'   observation) and `sd_y` (sd of the response over observations). The
#'   arguments `y` and `X` are then species means, and `se` is not used.
#' @return A named list.
#' @export
kernel_standata <- function(kernel, tree, y, X = NULL, se = NULL, use_prior = 1,
                          prior_sigma_scale = NULL, prior_beta_scale = NULL,
                          repeated = NULL) {
  ti <- tree_info(tree)
  L <- kernel_layout(kernel, ti)
  y <- as.numeric(y)
  if (length(y) != ti$N) stop("`y` must have one value per tip.", call. = FALSE)
  X <- if (is.null(X)) matrix(1, ti$N, 1) else as.matrix(X)
  sdy <- stats::sd(y)
  d <- c(list(N = ti$N, P = ncol(X), X = unname(X), y = y,
              se = se %||% rep(0, ti$N)),
         stan_tree_data(ti),
         list(use_prior = use_prior,
              prior_beta_scale = prior_beta_scale %||% (10 * max(1, sdy, abs(mean(y))))))
  if (L$dense) d$B <- tree_B(ti)
  if ("tC" %in% L$needs) d$tC <- tree_C(ti)
  if ("Dtree" %in% L$needs) d$Dtree <- tree_D(ti)
  if (L$quad_kind == "gh") {
    gh <- gauss_hermite(30)
    d <- c(d, list(n_gh = 30L, gh_x = gh$x, gh_w = gh$w))
  }
  if (L$quad_kind == "lam") {
    cmp <- L$kernel$components[[L$quad]]
    lq <- lambda_nodes(cmp$df)
    d <- c(d, list(n_lq = length(lq$lam), lq_lam = lq$lam, lq_w = lq$w))
  }
  for (j in seq_along(L$kernel$components)) {
    cmp <- L$kernel$components[[j]]; id <- L$ids[j]
    bv <- cmp$base_var(ti)
    if (!is.null(cmp$scale)) bv <- bv * mean(slope_values(cmp, ti)^2)
    d[[paste0(id, "_sigma_scale")]] <- prior_sigma_scale %||% (2 * sdy / sqrt(bv))
    for (p in names(cmp$pars)) {
      b <- resolve_bounds(cmp$pars[[p]], ti, cmp)
      d[[paste0(id, "_", p, "_lb")]] <- b[1]
      d[[paste0(id, "_", p, "_ub")]] <- b[2]
    }
    if (!is.null(cmp$standata)) d <- c(d, cmp$standata(ti, id, cmp))
  }
  if (!is.null(repeated))
    d <- c(d, slope_obs_data(L, repeated), list(M = length(repeated$species),
                   obs_sp = as.integer(repeated$species), y_obs = as.numeric(repeated$y),
                   X_obs = matrix(repeated$X, ncol = ncol(X)),
                   se_obs = repeated$se %||% rep(0, length(repeated$species)),
                   within_include = rep(1, ti$N),
                   within_sigma_scale = prior_sigma_scale %||% (2 * repeated$sd_y)))
  d
}

# x of each observation (centred) for varying slopes with several rows per species
slope_obs_data <- function(L, repeated) {
  j <- which(L$forms == "slope")
  if (!length(j)) return(list())
  stats::setNames(lapply(j, function(k) {
    sc <- L$kernel$components[[k]]$scale
    if (is.null(sc$x_obs)) stop("k_slope(): no values of x for the observations.", call. = FALSE)
    sc$x_obs - sc$center_value
  }), paste0(L$ids[j], "_xs_obs"))
}

# Gauss-Hermite nodes and weights (Golub-Welsch) for integrals against exp(-x^2)
gauss_hermite <- function(n) {
  J <- matrix(0, n, n)
  b <- sqrt(seq_len(n - 1) / 2)
  J[cbind(1:(n - 1), 2:n)] <- b
  J[cbind(2:n, 1:(n - 1))] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}

# Quadrature for a variance multiplier lam ~ inverse-gamma(df/2, df/2):
# trapezoid rule on log(lam) between extreme quantiles (smooth, decaying
# integrand, so the rule converges quickly). Returns nodes and weights that
# sum to one.
lambda_nodes <- function(df, n = 100) {
  b <- df / 2
  t <- seq(log(1 / stats::qgamma(1 - 1e-14, b, b)), log(1 / stats::qgamma(1e-14, b, b)),
           length.out = n)
  lam <- exp(t)
  w <- stats::dgamma(1 / lam, b, b) / lam      # density of log(lam)
  list(lam = lam, w = w / sum(w))
}
