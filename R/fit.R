# Fitting: phylo_gp() and methods for its results.

#' Fit a phylogenetic Gaussian-process model
#'
#' Fits `y ~ N(X beta, K)` where `K` is the covariance implied by `kernel` on
#' `tree`. The Stan program is written by [kernel_stancode()], compiled once
#' (and cached), and run with `cmdstanr`.
#'
#' @param formula Model formula, e.g. `brain ~ body`.
#' @param data Data frame. Species are matched to tip labels through
#'   `species` or, if `species` is `NULL`, the row names. With `species`, a
#'   species may have several rows (repeated measurements); see Details.
#' @param tree A `phylo` tree with branch lengths.
#' @param kernel A kernel, e.g. `k_bm() + k_white()`.
#' @param method `"sample"` (Bayesian, HMC) or `"optimize"` (maximum
#'   likelihood: no priors, no Jacobian, several starts). Kernels with latent
#'   rates (`k_relaxed()`, `k_rate_time()`) need `"sample"`.
#' @param species Optional name of the column of `data` holding species
#'   names. Needed when species have several rows.
#' @param se Optional name of a column with known measurement-error standard
#'   deviations.
#' @param prior_sigma_scale,prior_beta_scale Prior scales; see
#'   [kernel_standata()].
#' @param chains,iter_warmup,iter_sampling,adapt_delta,seed Passed to
#'   `cmdstanr`'s `$sample()`.
#' @param n_starts Number of optimisation starts for `method = "optimize"`.
#' @param ... Further arguments to `$sample()` or `$optimize()`.
#' @details
#' **Several rows per species.** When some species have more than one row,
#' the model is
#' \deqn{y_{ij} = x_{ij}^\top\beta + f_i + \epsilon_{ij}, \quad
#'   \epsilon_{ij} \sim N(0, \sigma_w^2),}
#' where \eqn{f} is the Gaussian process given by `kernel` over species and
#' \eqn{\sigma_w} (`within_sigma`) is the residual standard deviation within
#' species. A `k_white()` component is then the non-phylogenetic species
#' effect (the `(1 | species)` term of a mixed model). The likelihood is
#' exact and keeps its linear cost: it factorises into the species means,
#' which follow the kernel with an extra tip variance \eqn{\sigma_w^2/n_i},
#' and the deviations from those means. Predictors may vary within species.
#' Leave-one-out leaves out whole species. A varying slope
#' (`k_slope(x) * k_bm()`, latent) then applies to each row's own `x`, and
#' leave-one-out integrates the left-out species' slope over all its rows.
#'
#' With `se` and several rows per species, `se` is the known measurement
#' error of each row: \eqn{\epsilon_{ij} \sim N(0, \sigma_w^2 + se_{ij}^2)}.
#' The species means are then weighted by \eqn{1/(\sigma_w^2 + se_{ij}^2)},
#' and the likelihood stays exact and linear in the number of rows.
#' @return An object of class `phylogp_fit`.
#' @examples
#' \dontrun{
#' tr <- ape::rcoal(50)
#' d <- data.frame(y = kernel_simulate(k_bm() + k_white(), tr)[, 1])
#' fit <- phylo_gp(y ~ 1, d, tr, k_bm() + k_white())
#' fit
#' loo(fit)
#' }
#' @references
#' Rasmussen C.E., Williams C.K.I. 2006. Gaussian processes for machine learning. Cambridge, MA: MIT Press.
#'
#' Carpenter B., Gelman A., Hoffman M.D., Lee D., Goodrich B., Betancourt M., Brubaker M., Guo J., Li P., Riddell A. 2017. Stan: A probabilistic programming language. Journal of Statistical Software. 76:1–32.
#'
#' Felsenstein J. 1973. Maximum-likelihood estimation of evolutionary trees from continuous characters. American Journal of Human Genetics. 25:471–492.
#'
#' Ho L. si T., Ané C. 2014. A linear-time algorithm for Gaussian and non-Gaussian trait evolution models. Systematic Biology. 63:397–408.
#'
#' Ives A.R., Midford P.E., Garland T. 2007. Within-species variation and measurement error in phylogenetic comparative methods. Systematic Biology. 56:252–270.
#'
#' Hadfield J.D., Nakagawa S. 2010. General quantitative genetic methods for comparative biology: Phylogenies, taxonomies and multi-trait models for continuous and categorical characters. Journal of Evolutionary Biology. 23:494–508.
#' @export
phylo_gp <- function(formula, data, tree, kernel,
                     method = c("sample", "optimize"),
                     species = NULL, se = NULL,
                     prior_sigma_scale = NULL, prior_beta_scale = NULL,
                     chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                     adapt_delta = 0.9, seed = 1, n_starts = 8, ...) {
  method <- match.arg(method)
  kernel <- as_kernel(kernel)
  check_cmdstanr()
  md <- match_data(formula, data, tree, species, c(se, slope_columns(kernel)), se = se)
  if (!is.null(md$repeated)) {
    sl <- Filter(function(cmp) !is.null(cmp$scale), kernel$components)
    if (length(sl) > 1 || any(vapply(sl, function(cmp) !isTRUE(cmp$latent), NA)))
      stop("With several rows per species, varying slopes are available as one ",
           "k_slope(x) * k term with a branch-wise k (e.g. k_bm()) and latent = TRUE.",
           call. = FALSE)
  }
  kernel <- resolve_slopes(kernel, md$data, md$repeated$rows, species)
  ti <- tree_info(md$tree)
  L <- kernel_layout(kernel, ti)
  kernel <- L$kernel                       # with the regimes of k_regime() filled in
  if (method == "optimize" && L$latent)
    stop("Kernels with latent values (k_relaxed(), k_rate_time(), latent slopes, ",
         "k_white(df = )) need method = \"sample\".", call. = FALSE)

  code <- kernel_stancode(kernel, ti$tree, repeated = !is.null(md$repeated))
  standata <- kernel_standata(kernel, ti, md$y, md$X, md$se,
                            use_prior = as.integer(method == "sample"),
                            prior_sigma_scale = prior_sigma_scale,
                            prior_beta_scale = prior_beta_scale,
                            repeated = md$repeated)
  mod <- compile_phylogp(code)

  if (method == "sample") {
    fit <- mod$sample(data = standata, chains = chains, parallel_chains = chains,
                      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                      adapt_delta = adapt_delta, seed = seed, refresh = 0,
                      show_messages = FALSE, ...)
  } else {
    fits <- lapply(seq_len(n_starts), function(s) tryCatch(suppressWarnings({
      f <- mod$optimize(data = standata, seed = seed + s, jacobian = FALSE,
                        algorithm = "lbfgs", iter = 10000, refresh = 0,
                        show_messages = FALSE, show_exceptions = FALSE, ...)
      f$lp()
      f
    }), error = function(e) NULL))
    fits <- Filter(Negate(is.null), fits)
    if (!length(fits)) stop("All optimisation starts failed.", call. = FALSE)
    lp <- vapply(fits, function(f) f$lp(), numeric(1))
    fit <- fits[[which.max(lp)]]
  }
  structure(list(fit = fit, method = method, kernel = kernel, layout = L,
                 code = code, standata = standata, tree = ti$tree, tree_full = tree,
                 formula = formula, terms = md$terms, xlevels = md$xlevels,
                 species = species, data = md$data, coef_names = colnames(md$X),
                 repeated = md$repeated,
                 sampler = list(chains = chains, iter_warmup = iter_warmup,
                                iter_sampling = iter_sampling,
                                adapt_delta = adapt_delta, seed = seed)),
            class = "phylogp_fit")
}

check_cmdstanr <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    stop("Fitting needs the 'cmdstanr' package and CmdStan:\n",
         "  install.packages(\"cmdstanr\", repos = c(\"https://stan-dev.r-universe.dev\", getOption(\"repos\")))\n",
         "  cmdstanr::install_cmdstan()", call. = FALSE)
  if (is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)))
    stop("CmdStan not found; run cmdstanr::install_cmdstan().", call. = FALSE)
}

# Compile a generated program once; identical code reuses the cached binary.
compile_phylogp <- function(code) {
  dir <- getOption("phylogp.cache_dir", tools::R_user_dir("phylogp", "cache"))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(fileext = ".stan")
  writeLines(code, tmp)
  hash <- unname(tools::md5sum(tmp))
  file <- file.path(dir, paste0("phylogp_", substr(hash, 1, 12), ".stan"))
  if (!file.exists(file)) file.copy(tmp, file)
  cmdstanr::cmdstan_model(file, quiet = TRUE)
}

# Match data rows to tips; drop species missing from either, or with NAs.
# With several rows per species (only through `species`), y and X become
# species means and `repeated` holds what the within-species part of the
# likelihood needs: rows per species, the species of each row, and the
# response, design matrix and measurement error (`se`) of each row.
match_data <- function(formula, data, tree, species, extra_vars = NULL, se = NULL) {
  data <- as.data.frame(data)
  sp <- if (is.null(species)) rownames(data) else as.character(data[[species]])
  if (is.null(sp)) stop("Species names (row names or `species`) must be present.", call. = FALSE)
  if (anyDuplicated(sp) && is.null(species))
    stop("Row names must be unique; for several rows per species give the ",
         "species column in `species`.", call. = FALSE)
  vars <- unique(c(all.vars(formula), extra_vars))
  ok <- stats::complete.cases(data[, vars, drop = FALSE]) & !is.na(sp)
  keep <- intersect(tree$tip.label, unique(sp[ok]))
  n_drop_tree <- length(setdiff(tree$tip.label, keep))
  n_drop_rows <- sum(!(ok & sp %in% keep))
  if (length(keep) < 3)
    stop("Fewer than 3 species match the tree.",
         if (is.null(species)) " Species are taken from the row names; give the column of species names in `species`.",
         call. = FALSE)
  if (n_drop_tree) {
    message(sprintf("Dropping %d tips without data.", n_drop_tree))
    tree <- drop_tips(tree, setdiff(tree$tip.label, keep))
  }
  if (n_drop_rows) message(sprintf("Dropping %d rows not in the tree or with missing values.", n_drop_rows))
  ti <- tree_info(tree)
  rows <- data[ok & sp %in% keep, , drop = FALSE]
  sp_rows <- sp[ok & sp %in% keep]
  mf <- stats::model.frame(formula, rows)
  terms <- attr(mf, "terms")
  xlevels <- stats::.getXlevels(terms, mf)
  y_rows <- stats::model.response(mf)
  X_rows <- stats::model.matrix(formula, mf)
  if (!anyDuplicated(sp_rows)) {
    o <- match(ti$labels, sp_rows)
    data_sp <- rows[o, , drop = FALSE]
    rownames(data_sp) <- ti$labels
    return(list(tree = ti$tree, data = data_sp, y = unname(y_rows[o]),
                X = X_rows[o, , drop = FALSE], terms = terms, xlevels = xlevels,
                se = if (is.null(se)) NULL else data_sp[[se]], repeated = NULL))
  }
  idx <- match(sp_rows, ti$labels)
  n <- tabulate(idx, ti$N)
  y_mean <- as.numeric(rowsum(y_rows, idx, reorder = TRUE)) / n
  X_mean <- rowsum(X_rows, idx, reorder = TRUE) / n
  dimnames(X_mean) <- list(ti$labels, colnames(X_rows))
  data_sp <- rows[match(ti$labels, sp_rows), , drop = FALSE]
  rownames(data_sp) <- ti$labels
  list(tree = ti$tree, data = data_sp, y = y_mean, X = X_mean,
       terms = terms, xlevels = xlevels, se = NULL,
       repeated = list(n = n, species = idx, y = as.numeric(y_rows), X = unname(X_rows),
                       se = if (is.null(se)) NULL else as.numeric(rows[[se]]),
                       sd_y = stats::sd(y_rows), rows = rows))
}

# Names of the parameters that summarise a fit (not latent variables)
fit_par_names <- function(x) {
  L <- x$layout
  b <- paste0("coefs[", seq_along(x$coef_names), "]")
  cp <- unlist(lapply(seq_along(L$ids), function(j) {
    cmp <- L$kernel$components[[j]]
    id <- L$ids[j]
    c(paste0(id, "_sigma"), if (length(cmp$pars)) paste0(id, "_", names(cmp$pars)),
      if (length(cmp$hyper)) paste0(id, "_", cmp$hyper))
  }))
  c(b, cp, if (!is.null(x$repeated)) "within_sigma")
}

#' @export
print.phylogp_fit <- function(x, digits = 3, ...) {
  cat("Phylogenetic Gaussian process\n")
  cat("  Formula:", deparse(x$formula), "\n")
  cat("  Kernel: ", paste(vapply(x$kernel$components, `[[`, "", "label"), collapse = " + "), "\n")
  cat("  Tips:   ", length(x$tree$tip.label),
      if (!is.null(x$repeated)) sprintf("(%d rows; within-species sd: within_sigma)", length(x$repeated$species)), "\n")
  cat("  Likelihood:", if (x$layout$dense) "dense" else "pruning (linear time)", "\n")
  cat("  Method: ", if (x$method == "sample") "Bayesian (HMC)" else "maximum likelihood", "\n\n")
  s <- summary(x)
  print(s, digits = digits)
  if (x$method == "optimize") cat("\nlog-likelihood:", format(x$fit$lp(), digits = 10), "\n")
  invisible(x)
}

#' @export
summary.phylogp_fit <- function(object, ...) {
  pn <- fit_par_names(object)
  lab <- pn
  lab[seq_along(object$coef_names)] <- object$coef_names
  if (object$method == "optimize") {
    est <- object$fit$mle()[pn]
    out <- data.frame(parameter = lab, estimate = unname(est))
  } else {
    s <- as.data.frame(object$fit$summary(variables = pn, "mean", "sd",
                                          ~stats::quantile(.x, c(0.05, 0.95)),
                                          "rhat", "ess_bulk"))
    names(s) <- c("variable", "mean", "sd", "q5", "q95", "rhat", "ess_bulk")
    out <- data.frame(parameter = lab, s[, -1])
  }
  out
}

#' @export
coef.phylogp_fit <- function(object, ...) {
  s <- summary(object)
  stats::setNames(s[[2]], s$parameter)
}

#' @export
logLik.phylogp_fit <- function(object, ...) {
  if (object$method != "optimize")
    stop("logLik() is available for maximum-likelihood fits; use loo() for Bayesian fits.", call. = FALSE)
  structure(object$fit$lp(), df = length(fit_par_names(object)),
            nobs = length(object$tree$tip.label), class = "logLik")
}

#' Leave-one-out cross-validation for a phylogenetic GP fit
#'
#' For Bayesian fits, Pareto-smoothed importance-sampling LOO (`loo::loo()`)
#' built on the exact conditional densities `p(y_i | y_-i, theta)`, which the
#' generated Stan program computes in linear time for branch-wise kernels.
#' For maximum-likelihood fits, the conditional densities at the estimates.
#' Check the Pareto-k diagnostics: values above 0.7 flag influential species
#' whose contribution should be recomputed by refitting.
#'
#' For a fit across a sample of trees ([phylo_gp_trees()]), PSIS-LOO is run on
#' each tree and the left-out density of each species is averaged over trees
#' with equal weights, \eqn{\log \frac{1}{K}\sum_T p_T(y_i \mid y_{-i})}. The
#' per-tree results are kept in `attr(, "per_tree")`, and each species' Pareto
#' k is its largest value across trees.
#'
#' @param x A `phylogp_fit`, or a `phylogp_trees` object.
#' @param ... Passed to `loo::loo()`.
#' @return A `loo` object (Bayesian fits) or a vector of pointwise log
#'   densities (maximum-likelihood fits).
#' @importFrom loo loo
#' @references
#' Sundararajan S., Keerthi S.S. 2001. Predictive approaches for choosing hyperparameters in Gaussian processes. Neural Computation. 13:1103–1118.
#'
#' Vehtari A., Gelman A., Gabry J. 2017. Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing. 27:1413–1432.
#'
#' Bürkner P.-C., Gabry J., Vehtari A. 2021. Efficient leave-one-out cross-validation for Bayesian non-factorized normal and Student-t models. Computational Statistics. 36:1243–1261.
#' @exportS3Method
loo.phylogp_fit <- function(x, ...) {
  if (x$method == "optimize") {
    return(stats::setNames(as.numeric(x$fit$mle()[grep("^log_lik", names(x$fit$mle()))]),
                           x$tree$tip.label))
  }
  ll <- x$fit$draws("log_lik", format = "draws_array")
  r_eff <- loo::relative_eff(exp(ll))
  loo::loo(ll, r_eff = r_eff, ...)
}

#' @export
loo::loo

