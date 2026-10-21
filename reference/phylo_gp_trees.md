# Fit a model across a sample of trees

Fits the same model to each tree of a sample (for example, trees drawn
from the posterior of a dated phylogeny) and pools the results, so that
estimates and model comparisons integrate over phylogenetic uncertainty.
Each tree gets the same number of posterior draws, i.e. trees are
weighted equally, as in the common practice of repeating an analysis
over a tree sample. The Stan program is compiled once and reused for
every tree.

## Usage

``` r
phylo_gp_trees(
  formula,
  data,
  trees,
  kernel,
  n_trees = NULL,
  method = c("sample", "optimize"),
  species = NULL,
  se = NULL,
  seed = 1,
  ...
)
```

## Arguments

- formula, data, kernel, species, se:

  As in
  [`phylo_gp()`](https://thomaspueschel.github.io/phylogp/reference/phylo_gp.md).

- trees:

  A `multiPhylo` object or a list of `phylo` trees.

- n_trees:

  Optional number of trees to use, drawn at random from `trees`; all
  trees by default.

- method:

  `"sample"` (pooled posterior) or `"optimize"` (the distribution of
  maximum-likelihood estimates across trees).

- seed:

  Random seed for tree selection and fitting.

- ...:

  Further arguments to
  [`phylo_gp()`](https://thomaspueschel.github.io/phylogp/reference/phylo_gp.md),
  e.g. `chains`, `iter_sampling`.

## Value

An object of class `phylogp_trees` holding the individual fits.
[`summary()`](https://rdrr.io/r/base/summary.html) gives pooled
estimates and the spread of per-tree estimates;
[`loo()`](https://mc-stan.org/loo/reference/loo.html) gives
leave-one-out cross-validation for the pooled posterior;
[`coef()`](https://rdrr.io/r/stats/coef.html) gives per-tree estimates.

## References

Villemereuil P. de, Wells J.A., Edwards R.D., Blomberg S.P. 2012.
Bayesian models for comparative analysis integrating phylogenetic
uncertainty. BMC Evolutionary Biology. 12:102.
