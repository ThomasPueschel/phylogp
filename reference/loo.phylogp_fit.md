# Leave-one-out cross-validation for a phylogenetic GP fit

For Bayesian fits, Pareto-smoothed importance-sampling LOO
([`loo::loo()`](https://mc-stan.org/loo/reference/loo.html)) built on
the exact conditional densities `p(y_i | y_-i, theta)`, which the
generated Stan program computes in linear time for branch-wise kernels.
For maximum-likelihood fits, the conditional densities at the estimates.
Check the Pareto-k diagnostics: values above 0.7 flag influential
species whose contribution should be recomputed by refitting.

## Usage

``` r
# S3 method for class 'phylogp_fit'
loo(x, ...)

# S3 method for class 'phylogp_trees'
loo(x, ...)
```

## Arguments

- x:

  A `phylogp_fit`, or a `phylogp_trees` object.

- ...:

  Passed to [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html).

## Value

A `loo` object (Bayesian fits) or a vector of pointwise log densities
(maximum-likelihood fits).

## Details

For a fit across a sample of trees
([`phylo_gp_trees()`](https://thomaspueschel.github.io/phylogp/reference/phylo_gp_trees.md)),
PSIS-LOO is run on each tree and the left-out density of each species is
averaged over trees with equal weights, \\\log \frac{1}{K}\sum_T p_T(y_i
\mid y\_{-i})\\. The per-tree results are kept in `attr(, "per_tree")`,
and each species' Pareto k is its largest value across trees.

## References

Sundararajan S., Keerthi S.S. 2001. Predictive approaches for choosing
hyperparameters in Gaussian processes. Neural Computation. 13:1103–1118.

Vehtari A., Gelman A., Gabry J. 2017. Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. Statistics and
Computing. 27:1413–1432.

Bürkner P.-C., Gabry J., Vehtari A. 2021. Efficient leave-one-out
cross-validation for Bayesian non-factorized normal and Student-t
models. Computational Statistics. 36:1243–1261.
