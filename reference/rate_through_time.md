# Rate of evolution through time

Posterior summaries of the rate of evolution in each time bin of a
[`k_rate_time()`](https://thomaspueschel.github.io/phylogp/reference/k_rate_time.md)
component.

## Usage

``` r
rate_through_time(fit, probs = c(0.05, 0.95))
```

## Arguments

- fit:

  A Bayesian `phylogp_fit` whose kernel contains
  [`k_rate_time()`](https://thomaspueschel.github.io/phylogp/reference/k_rate_time.md).

- probs:

  Quantiles to report.

## Value

A data frame with one row per bin: start and end time (from the root),
posterior median and quantiles of the rate `sigma^2(t)`.
