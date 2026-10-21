# Ancestral states

Posterior distributions of the trait at the internal nodes of the tree
(kriging at internal points of the Gaussian process), computed in linear
time and averaged over posterior draws. Available for intercept-only
models whose kernel is a process along the whole tree: sums of
branch-wise components (plus residual or measurement error at the tips),
or a single Ornstein-Uhlenbeck component. For Pagel's kappa and delta,
reconstructions refer to the transformed tree; for Pagel's lambda, to
the tree with internal branches scaled by lambda.

## Usage

``` r
ancestral(object, ndraws = 500, probs = c(0.05, 0.95), tips = FALSE, seed = 1)
```

## Arguments

- object:

  A `phylogp_fit` of an intercept-only model (`y ~ 1`).

- ndraws:

  Number of posterior draws used.

- probs:

  Quantiles to report.

- tips:

  Also return the tips?

- seed:

  Random seed.

## Value

A data frame with one row per node (numbered as in `object$tree`):
posterior mean, standard deviation and quantiles.

## Details

For OU, reconstructions use the OU process itself, not its branch-wise
tip kernel: on an ultrametric tree OU and an accelerating early burst
have the same tip distribution but different ancestral states.

## References

Schluter D., Price T., Mooers A.Ø., Ludwig D. 1997. Likelihood of
ancestor states in adaptive radiation. Evolution. 51:1699–1711.

Garland T. Jr., Ives A.R. 2000. Using the past to predict the present:
Confidence intervals for regression equations in phylogenetic
comparative methods. The American Naturalist. 155:346–364.
