# Plot branch rates on the tree

Colours every branch by its posterior median relative rate (log scale)
from a relaxed-rate fit: blue slower, red faster than a typical branch.

## Usage

``` r
plot_edge_rates(fit, ...)
```

## Arguments

- fit:

  A Bayesian `phylogp_fit` whose kernel contains
  [`k_relaxed()`](https://thomaspueschel.github.io/phylogp/reference/k_relaxed.md).

- ...:

  Further arguments to
  [`plot_tree_rates()`](https://thomaspueschel.github.io/phylogp/reference/plot_tree_rates.md).

## Value

The output of
[`edge_rates()`](https://thomaspueschel.github.io/phylogp/reference/edge_rates.md),
invisibly.
