# Species- and node-level slopes from a varying-slope fit

Posterior summaries of the regression slope at every node of the tree,
`beta + v`, from a model with a latent `k_slope(x) * k` term. The mean
slope `beta` is the coefficient of `x` in the model formula.

## Usage

``` r
varying_slopes(fit, probs = c(0.05, 0.95))
```

## Arguments

- fit:

  A Bayesian `phylogp_fit` with a `k_slope(x) * ...` term fitted with
  `latent = TRUE`.

- probs:

  Quantiles to report.

## Value

A data frame with one row per node of `fit$tree` (tips first): node
number, tip label (`NA` for internal nodes), posterior mean, standard
deviation and quantiles of the slope.
