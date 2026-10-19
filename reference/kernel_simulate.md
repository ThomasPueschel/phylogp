# Simulate traits from a kernel

Simulate traits from a kernel

## Usage

``` r
kernel_simulate(kernel, tree, pars = NULL, nsim = 1, mean = 0)
```

## Arguments

- kernel:

  A kernel, e.g. `k_bm() + k_white()`.

- tree:

  A `phylo` tree (or the output of
  [`tree_info()`](https://thomaspueschel.github.io/phylogp/reference/tree_info.md)).

- pars:

  Optional parameter values, as a list named by component id (see
  `print(kernel)`), e.g.
  `list(lambda1 = list(sigma = 1, lambda = 0.3))`. Unspecified values
  take the defaults of each component. For components with latent rates
  ([`k_relaxed()`](https://thomaspueschel.github.io/phylogp/reference/k_relaxed.md),
  [`k_rate_time()`](https://thomaspueschel.github.io/phylogp/reference/k_rate_time.md)),
  a rate path is drawn from the prior unless supplied (`z` or
  `lograte`).

- nsim:

  Number of simulated data sets.

- mean:

  Mean of the trait (a number or a vector in tip order).

## Value

An N x nsim matrix with tip labels as row names.

## Examples

``` r
tr <- ape::rcoal(20)
y <- kernel_simulate(k_ou(), tr, pars = list(ou1 = list(alpha = 3)))
```
