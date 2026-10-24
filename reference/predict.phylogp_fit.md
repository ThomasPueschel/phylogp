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
  [`k_white()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md)
  and, for fits with several rows per species, the within-species
  variation) or `"latent"` (the evolutionary process at the tip).

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

## Details

For kernels with latent values, the values a new species needs are drawn
for each posterior draw from their conditional prior given the fitted
ones: relaxed log rates and latent slope deviations by a Brownian bridge
along the fitted branch the species joins and Brownian motion along its
own branch, a Student-t variance multiplier from its inverse-gamma
prior; rates through time keep the time bins of the fitted tree.

## References

Garland T. Jr., Ives A.R. 2000. Using the past to predict the present:
Confidence intervals for regression equations in phylogenetic
comparative methods. The American Naturalist. 155:346–364.

Rasmussen C.E., Williams C.K.I. 2006. Gaussian processes for machine
learning. Cambridge, MA: MIT Press.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- phylo_gp(y ~ x, d, tree, k_bm() + k_white())
predict(fit, newdata = data.frame(x = 1.2, row.names = "missing_species"))
} # }
```
