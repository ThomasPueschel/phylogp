# Fitting: phylo_gp() and methods for its results.

#' Fit a phylogenetic Gaussian-process model
#'
#' Fits `y ~ N(X beta, K)` where `K` is the covariance implied by `kernel` on
#' `tree`. The Stan program is written by [kernel_stancode()], compiled once
#' (and cached), and run with `cmdstanr`.
#'
#' @param formula Model formula, e.g. `brain ~ body`.
#' @param data Data frame. Species are matched to tip labels through
#'   `species` or, if `species` is `NULL`, the row names.
#' @param tree A `phylo` tree with branch lengths.
#' @param kernel A kernel, e.g. `k_bm() + k_white()`.
#' @param method `"sample"` (Bayesian, HMC) or `"optimize"` (maximum
#'   likelihood: no priors, no Jacobian, several starts). Kernels with latent
#'   rates (`k_relaxed()`, `k_rate_time()`) need `"sample"`.
#' @param species Optional name of the column of `data` holding species names.
#' @param se Optional name of a column with known measurement-error standard
#'   deviations.
#' @param prior_sigma_scale,prior_beta_scale Prior scales; see
#'   [kernel_standata()].
#' @param chains,iter_warmup,iter_sampling,adapt_delta,seed Passed to
#'   `cmdstanr`'s `$sample()`.
#' @param n_starts Number of optimisation starts for `method = "optimize"`.
#' @param ... Further arguments to `$sample()` or `$optimize()`.
#' @return An object of class `phylogp_fit`.
#' @examples
#' \dontrun{
#' tr <- ape::rcoal(50)
#' d <- data.frame(y = kernel_simulate(k_bm() + k_white(), tr)[, 1])
#' fit <- phylo_gp(y ~ 1, d, tr, k_bm() + k_white())
#' fit
#' loo(fit)
#' }
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
  kernel <- resolve_slopes(kernel, md$data)
  ti <- tree_info(md$tree)
  L <- kernel_layout(kernel, ti)
  if (method == "optimize" && L$latent)
    stop("Kernels with latent values (k_relaxed(), k_rate_time(), latent slopes, ",
         "k_white(df = )) need method = \"sample\".", call. = FALSE)

  code <- kernel_stancode(kernel, ti$tree)
  standata <- kernel_standata(kernel, ti, md$y, md$X, md$se,
                            use_prior = as.integer(method == "sample"),
                            prior_sigma_scale = prior_sigma_scale,
                            prior_beta_scale = prior_beta_scale)
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

# Match data rows to tips; drop species missing from either, or with NAs
match_data <- function(formula, data, tree, species, extra_vars = NULL, se = NULL) {
  data <- as.data.frame(data)
  sp <- if (is.null(species)) rownames(data) else as.character(data[[species]])
  if (is.null(sp) || anyDuplicated(sp))
    stop("Species names (row names or `species`) must be present and unique.", call. = FALSE)
  rownames(data) <- sp
  vars <- unique(c(all.vars(formula), extra_vars))
  ok <- stats::complete.cases(data[, vars, drop = FALSE])
  keep <- intersect(tree$tip.label, sp[ok])
  n_drop_tree <- length(setdiff(tree$tip.label, keep))
  n_drop_data <- length(setdiff(sp, keep))
  if (length(keep) < 3) stop("Fewer than 3 species match the tree.", call. = FALSE)
  if (n_drop_tree) {
    message(sprintf("Dropping %d tips without data.", n_drop_tree))
    tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
  }
  if (n_drop_data) message(sprintf("Dropping %d rows not in the tree or with missing values.", n_drop_data))
  ti <- tree_info(tree)
  data <- data[ti$labels, , drop = FALSE]
  mf <- stats::model.frame(formula, data)
  list(tree = ti$tree, data = data, y = stats::model.response(mf),
       X = stats::model.matrix(formula, mf),
       terms = attr(mf, "terms"), xlevels = stats::.getXlevels(attr(mf, "terms"), mf),
       se = if (is.null(se)) NULL else data[[se]])
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
  c(b, cp)
}

#' @export
print.phylogp_fit <- function(x, digits = 3, ...) {
  cat("Phylogenetic Gaussian process\n")
  cat("  Formula:", deparse(x$formula), "\n")
  cat("  Kernel: ", paste(vapply(x$kernel$components, `[[`, "", "label"), collapse = " + "), "\n")
  cat("  Tips:   ", length(x$tree$tip.label), "\n")
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
#' @param x A `phylogp_fit`.
#' @param ... Passed to `loo::loo()`.
#' @return A `loo` object (Bayesian fits) or a vector of pointwise log
#'   densities (maximum-likelihood fits).
#' @importFrom loo loo
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

