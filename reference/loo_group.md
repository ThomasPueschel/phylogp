# Leave-clade-out cross-validation

Leave-one-out mostly asks how well a species is predicted from its
closest relatives. Leaving out whole clades (or any other groups, such
as orders or families) asks a harder question that bears more directly
on the evolutionary process: how well does the model predict a group of
species from the rest of the tree?

## Usage

``` r
loo_group(fit, groups, ndraws = 400, ...)
```

## Arguments

- fit:

  A Bayesian `phylogp_fit`.

- groups:

  A vector giving the group of each species, named by species (tip
  label), or an unnamed vector in the order of the data rows used in the
  fit (`rownames(fit$data)`).

- ndraws:

  Number of posterior draws of each refit used.

- ...:

  Further arguments to `$sample()`.

## Value

An object of class `c("kfold", "loo")` usable with
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html):
estimates of the expected log predictive density summed over groups
(`elpd_kfold`), and a `pointwise` matrix with one row per group. The
attribute `groups` gives group names and sizes, and `species_elpd` the
marginal predictive density of each left-out species.

## Details

For each group the model is refitted without the group's data (its
species are given a very large measurement error, so the tree and any
latent values on it stay the same), and the joint predictive density of
the group, `log mean_s p(y_group | y_rest, theta_s)`, is computed. For
branch-wise kernels this is the difference of two pruning likelihoods,
so it costs O(N) per draw; dense kernels use dense conditioning.

## References

Vehtari A., Gelman A., Gabry J. 2017. Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. Statistics and
Computing. 27:1413–1432.

## Examples

``` r
if (FALSE) { # \dontrun{
lco <- loo_group(fit, groups = setNames(d$order, rownames(d)))
loo::loo_compare(lco_model_a, lco_model_b)
} # }
```
