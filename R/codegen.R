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

#' Stan code for a kernel
#'
#' Writes the complete Stan program that [phylo_gp()] compiles for a kernel.
#' The program can be inspected, edited and run directly with `cmdstanr`
#' together with the data from [make_standata()].
#'
#' @param kernel A kernel, e.g. `k_bm() + k_white()`.
#' @param tree Optional tree. It decides whether components such as
#'   `k_ou()` can use the branch-wise form (ultrametric trees) or need dense
#'   algebra. Without a tree an ultrametric tree is assumed.
#' @return A character string with class `phylogp_stancode`.
#' @examples
#' cat(make_stancode(k_bm() + k_white()))
#' @export
make_stancode <- function(kernel, tree = NULL) {
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
    if (L$quad_kind == "lam") c("int<lower=1> n_lq;", "vector<lower=0>[n_lq] lq_lam;", "vector<lower=0>[n_lq] lq_w;"))
  for (j in seq_along(comps)) {
    cmp <- comps[[j]]; id <- L$ids[j]
    data <- c(data, sprintf("// %s: %s", id, cmp$label),
              sprintf("real<lower=0> %s_sigma_scale;", id))
    for (p in names(cmp$pars))
      data <- c(data, sprintf("real %1$s_%2$s_lb;", id, p),
                sprintf("real<lower=%1$s_%2$s_lb> %1$s_%2$s_ub;", id, p))
    data <- c(data, sub_id(cmp$stan_data, id))
  }

  # ---- parameters ----
  pars <- "vector[P] coefs;              // regression coefficients"
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
  asm <- c("vector[E] w = rep_vector(0, E);", "real rv = 0;", "vector[N] t2 = square(se);")
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
  asm <- c(asm, "vector[N] r = y - X * coefs;",
           unlist(lapply(seq_along(comps), function(j) sub_id(comps[[j]]$stan_shift, L$ids[j]))))

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
  model <- c("{", ind(c(asm, lik)), "}", latent_prior,
             "if (use_prior == 1) {", ind(priors), "}")

  # ---- generated quantities: leave-one-out densities p(y_i | y_-i, theta) ----
  # mo: mean and variance of r_i given the other species. A left-out species'
  # own latent values (its slope or rate on the tip edge) are informed only by
  # y_i, so they are integrated out against their conditional prior: in closed
  # form (latent slopes) or by quadrature (relaxed rates, heavy-tailed residuals).
  loo <- c(if (L$dense) "matrix[N, 2] mo = gaussian_loo_moments(r, K);" else
             "matrix[N, 2] mo = tree_loo_moments(r, w, rv, t2, anc, des, n_nodes, root);",
           "for (i in 1:N) {")
  body <- c("int e = tip_edge[i];", "real ri = r[i];",
            if (L$dense) "real vi = mo[i, 2];" else "real vi = mo[i, 2] + t2[i];")
  for (j in L$closed) {
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
  loo <- c(loo, ind(body), "}")
  gq <- c("vector[N] log_lik;", "{", ind(c(asm, loo)), "}")

  fun_file <- system.file("stan", "phylogp_functions.stan", package = "phylogp")
  funs <- readLines(fun_file, warn = FALSE)
  header <- c(
    "// Generated by the phylogp R package (phylogenetic Gaussian processes)",
    sprintf("// Kernel: %s", paste(vapply(comps, `[[`, "", "label"), collapse = " + ")),
    sprintf("// Likelihood: %s", if (L$dense) "dense (O(N^3))" else "pruning (O(N))"))
  code <- c(header,
            "functions {", funs, "}",
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

#' Data for the Stan program of a kernel
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
#' @return A named list.
#' @export
make_standata <- function(kernel, tree, y, X = NULL, se = NULL, use_prior = 1,
                          prior_sigma_scale = NULL, prior_beta_scale = NULL) {
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
  d
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
