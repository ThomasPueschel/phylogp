# Predict trait values for species

Posterior predictive distributions of the trait for species in the tree:
species that were fitted (smoothed values) or species in the tree
without trait data (phylogenetic imputation, i.e. kriging on the tree).
Each prediction conditions on all observed species, through the kernel
and the regression part of the model, and integrates over the posterior
draws.

## Usage

``` r
# S3 method for class 'phylogp_fit'
predict(
  object,
  newdata = NULL,
  tree = NULL,
  type = c("response", "latent"),
  ndraws = 500,
  probs = c(0.05, 0.95),
  seed = 1,
  ...
)
```

## Arguments

- object:

  A `phylogp_fit`.

- newdata:

  Optional data frame with the predictor values of the species to
  predict, named as in the fit (row names or the `species` column). Its
  species must be in `tree`. If `NULL`, the fitted species are
  predicted.

- tree:

  Tree containing the species to predict; by default the tree given to
  [`phylo_gp()`](https://thomaspueschel.github.io/phylogp/reference/phylo_gp.md).

- type:

  `"response"` (a new observation of the trait, including residual
  variation from
  [`k_white()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md))
  or `"latent"` (the evolutionary process at the tip).

- ndraws:

  Number of posterior draws used.

- probs:

  Quantiles to report.

- seed:

  Random seed for draw selection and quantiles.

- ...:

  Not used.

## Value

A data frame with one row per species: whether it was observed,
posterior mean, standard deviation and quantiles.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- phylo_gp(y ~ x, d, tree, k_bm() + k_white())
predict(fit, newdata = data.frame(x = 1.2, row.names = "missing_species"))
} # }
```
