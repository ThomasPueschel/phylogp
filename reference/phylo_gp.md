# Fit a phylogenetic Gaussian-process model

Fits `y ~ N(X beta, K)` where `K` is the covariance implied by `kernel`
on `tree`. The Stan program is written by
[`kernel_stancode()`](https://thomaspueschel.github.io/phylogp/reference/kernel_stancode.md),
compiled once (and cached), and run with `cmdstanr`.

## Usage

``` r
phylo_gp(
  formula,
  data,
  tree,
  kernel,
  method = c("sample", "optimize"),
  species = NULL,
  se = NULL,
  prior_sigma_scale = NULL,
  prior_beta_scale = NULL,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.9,
  seed = 1,
  n_starts = 8,
  ...
)
```

## Arguments

- formula:

  Model formula, e.g. `brain ~ body`.

- data:

  Data frame. Species are matched to tip labels through `species` or, if
  `species` is `NULL`, the row names.

- tree:

  A `phylo` tree with branch lengths.

- kernel:

  A kernel, e.g. `k_bm() + k_white()`.

- method:

  `"sample"` (Bayesian, HMC) or `"optimize"` (maximum likelihood: no
  priors, no Jacobian, several starts). Kernels with latent rates
  ([`k_relaxed()`](https://thomaspueschel.github.io/phylogp/reference/k_relaxed.md),
  [`k_rate_time()`](https://thomaspueschel.github.io/phylogp/reference/k_rate_time.md))
  need `"sample"`.

- species:

  Optional name of the column of `data` holding species names.

- se:

  Optional name of a column with known measurement-error standard
  deviations.

- prior_sigma_scale, prior_beta_scale:

  Prior scales; see
  [`kernel_standata()`](https://thomaspueschel.github.io/phylogp/reference/kernel_standata.md).

- chains, iter_warmup, iter_sampling, adapt_delta, seed:

  Passed to `cmdstanr`'s `$sample()`.

- n_starts:

  Number of optimisation starts for `method = "optimize"`.

- ...:

  Further arguments to `$sample()` or `$optimize()`.

## Value

An object of class `phylogp_fit`.

## Examples

``` r
if (FALSE) { # \dontrun{
tr <- ape::rcoal(50)
d <- data.frame(y = kernel_simulate(k_bm() + k_white(), tr)[, 1])
fit <- phylo_gp(y ~ 1, d, tr, k_bm() + k_white())
fit
loo(fit)
} # }
```
