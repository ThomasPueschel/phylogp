# Rates that depend on mapped regimes

Brownian motion whose rate depends on a discrete regime painted on the
tree, as in multi-rate Brownian motion models (O'Meara et al. 2006). The
regimes come from the tree: give
[`phylo_gp()`](https://thomaspueschel.github.io/phylogp/reference/phylo_gp.md)
a stochastic character map (a `simmap` from
[`phytools::make.simmap()`](https://rdrr.io/pkg/phytools/man/make.simmap.html),
or a tree painted with
[`phytools::paintSubTree()`](https://rdrr.io/pkg/phytools/man/paintSubTree.html)),
whose `mapped.edge` holds the time each branch spends in each regime.
Regimes may change along a branch.

## Usage

``` r
k_regime(
  states = NULL,
  reference = NULL,
  rate_scale = 1,
  lower = 0.001,
  upper = 1000
)
```

## Arguments

- states:

  Optional regime names. By default, all regimes mapped on the tree
  (sorted).

- reference:

  Regime whose rate is 1; by default the first of `states`.

- rate_scale:

  Standard deviation of the lognormal priors of the relative rates.

- lower, upper:

  Bounds of the relative rates.

## Value

A kernel component.

## Details

Each edge accumulates variance `sigma^2 * sum_r rate_r * t_r`, where
`t_r` is the time it spends in regime `r`. The rate of the reference
regime is 1, so `sigma` is its standard deviation per unit time and the
`rate_` parameters are the rates of the other regimes relative to it.
The kernel is branch-wise, so it keeps the linear-time likelihood,
leave-one-out, prediction and ancestral states, and it can be fitted by
maximum likelihood. In Bayesian fits the relative rates get lognormal
priors centred on 1.

To integrate over uncertainty in the history of the regimes, fit a
sample of stochastic maps (a `multiSimmap`) with
[`phylo_gp_trees()`](https://thomaspueschel.github.io/phylogp/reference/phylo_gp_trees.md).

## References

O'Meara B.C., Ané C., Sanderson M.J., Wainwright P.C. 2006. Testing for
different rates of continuous trait evolution using likelihood.
Evolution. 60:922–933.

Huelsenbeck J.P., Nielsen R., Bollback J.P. 2003. Stochastic mapping of
morphological characters. Systematic Biology. 52:131–158.

Revell L.J. 2012. phytools: An R package for phylogenetic comparative
biology (and other things). Methods in Ecology and Evolution. 3:217–223.

## Examples

``` r
if (FALSE) { # \dontrun{
tr <- phytools::pbtree(n = 100, scale = 1)
Q <- matrix(c(-1, 1, 1, -1), 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
smap <- phytools::sim.history(tr, Q)
fit <- phylo_gp(y ~ 1, d, smap, k_regime() + k_white(), method = "optimize")
} # }
```
