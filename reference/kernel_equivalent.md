# Do two kernels imply the same distribution for the tips?

Different evolutionary processes can imply exactly, or nearly, the same
distribution for the species at the tips of a tree, and then no amount
of extant data can tell them apart. On an ultrametric tree, for example,
an Ornstein-Uhlenbeck process with a fixed root has the same kernel as
an early-burst model whose accelerating rate is twice the OU `alpha`.

## Usage

``` r
kernel_equivalent(
  k1,
  k2,
  tree,
  pars = NULL,
  n_grid = 8,
  mean = c("free", "zero"),
  n_starts = 5,
  seed = 1
)
```

## Arguments

- k1:

  Kernel whose tip distribution is to be reproduced.

- k2:

  Candidate kernel.

- tree:

  A `phylo` tree (or the output of
  [`tree_info()`](https://thomaspueschel.github.io/phylogp/reference/tree_info.md)).

- pars:

  Optional parameter values of `k1` (as in
  [`kernel_matrix()`](https://thomaspueschel.github.io/phylogp/reference/kernel_matrix.md));
  by default a grid of `n_grid` settings spanning the bounds of its
  parameters.

- n_grid:

  Number of parameter settings of `k1` tried.

- mean:

  `"free"` (unknown mean, contrasts only) or `"zero"`.

- n_starts:

  Number of starting points of the search over `k2`.

- seed:

  Random seed for the starting points.

## Value

A data frame (class `phylogp_equivalence`) with one row per setting of
`k1`: its parameters, the smallest KL found, the best parameters of
`k2`, and a verdict (`"identical"`, `"nearly identical"` or
`"distinguishable"`).

## Details

For each of several parameter settings of `k1`, `kernel_equivalent()`
searches over all parameters of `k2` (including the scales `sigma`) for
the tip distribution closest to that of `k1`, measured by the
Kullback-Leibler divergence \$\$KL = E\_{k1}\[\log p\_{k1}(y) - \log
p\_{k2}(y)\].\$\$ KL is the expected log-likelihood ratio in favour of
`k1` for one data set of these tips simulated from `k1`, when `k2`'s
parameters are chosen as well as possible. Values near 0 mean `k2` can
reproduce `k1`; values below about 1 mean the two will rarely be
distinguishable with data of this size.

With `mean = "free"` (default), the trait mean is an unknown parameter,
so only contrasts between species are compared, as in REML. Two kernels
that differ only by a constant added to every covariance (such as a
random versus a fixed root) are then equivalent. With `mean = "zero"`
the full covariance matrices are compared.

The search is numerical, with several starting points. A small KL is
therefore evidence of equivalence, but a large one means only that no
equivalent parameter values were found.

## Examples

``` r
tr <- ape::rcoal(30)
# OU with a fixed root vs accelerating early burst: identical on this tree
kernel_equivalent(k_ou(root = "fixed"), k_eb(upper = 20), tr, n_grid = 3)
#> Can eb() reproduce the tip distribution of ou(root = "fixed")?
#> 30 tips, ultrametric tree, unknown (contrasts compared) mean
#> 
#>  k1.ou1.sigma k1.ou1.alpha       KL k2.eb1.sigma k2.eb1.rate k2_at_bound
#>             1       0.0409 7.20e-10     9.51e-01      0.0819       FALSE
#>             1       1.2900 4.33e-09     2.06e-01      2.5900       FALSE
#>             1      40.9000 2.48e+00     3.27e-06     20.0000        TRUE
#>          verdict
#>        identical
#>        identical
#>  distinguishable
#> 
#> KL: expected log-likelihood ratio in favour of the first kernel for one
#> data set simulated from it, with the second kernel's parameters chosen as
#> well as possible.
#> In some settings the best parameters of the second kernel are at their bounds;
#> widen the bounds (e.g. k_eb(upper = )) before concluding that the kernels differ.
#> Equivalence depends on the parameter values of the first kernel. 
```
