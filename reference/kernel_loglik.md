# Log-likelihood of a kernel in R

Evaluates the Gaussian log-likelihood at given parameter values, with
the O(N) pruning recursion when the kernel is branch-wise and a Cholesky
factorisation otherwise. Useful for checking a kernel or a fit.

## Usage

``` r
kernel_loglik(kernel, tree, y, pars = NULL, mean = 0)
```

## Arguments

- kernel:

  A kernel, e.g. `k_bm() + k_white()`.

- tree:

  A `phylo` tree (or the output of
  [`tree_info()`](https://thomaspueschel.github.io/phylogp/reference/tree_info.md)).

- y:

  Trait values, named by tip label or in tip order.

- pars:

  Optional parameter values, as a list named by component id (see
  `print(kernel)`), e.g.
  `list(lambda1 = list(sigma = 1, lambda = 0.3))`. Unspecified values
  take the defaults of each component. For components with latent rates
  ([`k_relaxed()`](https://thomaspueschel.github.io/phylogp/reference/k_relaxed.md),
  [`k_rate_time()`](https://thomaspueschel.github.io/phylogp/reference/k_rate_time.md)),
  a rate path is drawn from the prior unless supplied (`z` or
  `lograte`).

- mean:

  Mean vector (or number) to subtract.

## Value

The log-likelihood.
