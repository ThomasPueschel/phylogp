# Exact leave-one-out for species with high Pareto k

PSIS-LOO approximates leaving a species out by reweighting the posterior
draws of the full fit. For influential species (Pareto k above about
0.7) the approximation is unreliable. `loo_refit()` refits the model
without each such species and replaces its contribution by the exact
value `log mean_s p(y_i | y_-i, theta_s)`, with `theta_s` drawn from the
posterior given the other species.

## Usage

``` r
loo_refit(fit, loo = NULL, k_threshold = 0.7, ...)
```

## Arguments

- fit:

  A Bayesian `phylogp_fit`.

- loo:

  Optional result of `loo(fit)`; computed if missing.

- k_threshold:

  Species with Pareto k above this value are refitted.

- ...:

  Further arguments to `$sample()`.

## Value

The `loo` object with the pointwise and total estimates updated for the
refitted species (their Pareto k set to 0, as in `loo::reloo()`), and an
attribute `refitted` naming them.

## Details

A species is left out by giving it a measurement-error standard
deviation many times larger than the spread of the data, so its value
carries no information while the tree, and any latent rates on it, stay
the same. Its predictive density is then computed in linear time
(branch-wise kernels) or by dense conditioning.
