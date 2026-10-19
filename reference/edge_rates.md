# Branch rates from a relaxed-rate fit

Posterior summaries of the rate of evolution on every edge, from a model
with a
[`k_relaxed()`](https://thomaspueschel.github.io/phylogp/reference/k_relaxed.md)
component. Rates are relative to the median edge rate of each posterior
draw, so 2 means twice as fast as a typical branch.

## Usage

``` r
edge_rates(fit, probs = c(0.05, 0.95))
```

## Arguments

- fit:

  A Bayesian `phylogp_fit` whose kernel contains
  [`k_relaxed()`](https://thomaspueschel.github.io/phylogp/reference/k_relaxed.md).

- probs:

  Quantiles to report.

## Value

A data frame with one row per edge, in the edge order of `fit$tree`,
giving the parent and child node, the posterior median and the requested
quantiles of the relative rate.
