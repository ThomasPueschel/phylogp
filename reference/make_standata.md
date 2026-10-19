# Data for the Stan program of a kernel

Data for the Stan program of a kernel

## Usage

``` r
make_standata(
  kernel,
  tree,
  y,
  X = NULL,
  se = NULL,
  use_prior = 1,
  prior_sigma_scale = NULL,
  prior_beta_scale = NULL
)
```

## Arguments

- kernel:

  A kernel.

- tree:

  A `phylo` tree, or the output of
  [`tree_info()`](https://thomaspueschel.github.io/phylogp/reference/tree_info.md).

- y:

  Response vector in the tip order of `tree`.

- X:

  Design matrix in the same order (default: intercept only).

- se:

  Known measurement-error standard deviations (default 0).

- use_prior:

  `0` for maximum likelihood (no priors), `1` for Bayesian.

- prior_sigma_scale:

  Optional scale of the half-normal priors on the component scales; by
  default `2 sd(y) / sqrt(base variance)`.

- prior_beta_scale:

  Optional scale of the normal priors on the regression coefficients.

## Value

A named list.
