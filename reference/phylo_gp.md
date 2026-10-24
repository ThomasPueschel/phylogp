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
  `species` is `NULL`, the row names. With `species`, a species may have
  several rows (repeated measurements); see Details.

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

  Optional name of the column of `data` holding species names. Needed
  when species have several rows.

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

## Details

**Several rows per species.** When some species have more than one row,
the model is \$\$y\_{ij} = x\_{ij}^\top\beta + f_i + \epsilon\_{ij},
\quad \epsilon\_{ij} \sim N(0, \sigma_w^2),\$\$ where \\f\\ is the
Gaussian process given by `kernel` over species and \\\sigma_w\\
(`within_sigma`) is the residual standard deviation within species. A
[`k_white()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md)
component is then the non-phylogenetic species effect (the
`(1 | species)` term of a mixed model). The likelihood is exact and
keeps its linear cost: it factorises into the species means, which
follow the kernel with an extra tip variance \\\sigma_w^2/n_i\\, and the
deviations from those means. Predictors may vary within species.
Leave-one-out leaves out whole species. A varying slope
(`k_slope(x) * k_bm()`, latent) then applies to each row's own `x`, and
leave-one-out integrates the left-out species' slope over all its rows.

With `se` and several rows per species, `se` is the known measurement
error of each row: \\\epsilon\_{ij} \sim N(0, \sigma_w^2 +
se\_{ij}^2)\\. The species means are then weighted by \\1/(\sigma_w^2 +
se\_{ij}^2)\\, and the likelihood stays exact and linear in the number
of rows.

## References

Rasmussen C.E., Williams C.K.I. 2006. Gaussian processes for machine
learning. Cambridge, MA: MIT Press.

Carpenter B., Gelman A., Hoffman M.D., Lee D., Goodrich B., Betancourt
M., Brubaker M., Guo J., Li P., Riddell A. 2017. Stan: A probabilistic
programming language. Journal of Statistical Software. 76:1–32.

Felsenstein J. 1973. Maximum-likelihood estimation of evolutionary trees
from continuous characters. American Journal of Human Genetics.
25:471–492.

Ho L. si T., Ané C. 2014. A linear-time algorithm for Gaussian and
non-Gaussian trait evolution models. Systematic Biology. 63:397–408.

Ives A.R., Midford P.E., Garland T. 2007. Within-species variation and
measurement error in phylogenetic comparative methods. Systematic
Biology. 56:252–270.

Hadfield J.D., Nakagawa S. 2010. General quantitative genetic methods
for comparative biology: Phylogenies, taxonomies and multi-trait models
for continuous and categorical characters. Journal of Evolutionary
Biology. 23:494–508.

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
