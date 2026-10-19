# Covariance matrix implied by a kernel

Covariance matrix implied by a kernel

## Usage

``` r
kernel_matrix(kernel, tree, pars = NULL)
```

## Arguments

- kernel:

  A kernel, e.g. `k_bm() + k_white()`.

- tree:

  A `phylo` tree (or the output of
  [`tree_info()`](https://thomaspueschel.github.io/phylogp/reference/tree_info.md)).

- pars:

  Optional parameter values, as a list named by component id (see
  `print(kernel)`), e.g.
  `list(lambda1 = list(sigma = 1, lambda = 0.3))`. Unspecified values
  take the defaults of each component. For components with latent rates
  ([`k_relaxed()`](https://thomaspueschel.github.io/phylogp/reference/k_relaxed.md),
  [`k_rate_time()`](https://thomaspueschel.github.io/phylogp/reference/k_rate_time.md)),
  a rate path is drawn from the prior unless supplied (`z` or
  `lograte`).

## Value

An N x N covariance matrix with tip labels as dimnames.

## Examples

``` r
tr <- ape::rcoal(10)
kernel_matrix(k_bm() + k_white(), tr)
#>           t4       t7       t5       t8      t10       t9       t6       t1
#> t4  2.470506 1.310641 1.090563 0.000000 0.000000 0.000000 0.000000 0.000000
#> t7  1.310641 2.470506 1.090563 0.000000 0.000000 0.000000 0.000000 0.000000
#> t5  1.090563 1.090563 2.470506 0.000000 0.000000 0.000000 0.000000 0.000000
#> t8  0.000000 0.000000 0.000000 2.470506 1.449474 1.449474 1.022726 1.022726
#> t10 0.000000 0.000000 0.000000 1.449474 2.470506 1.460181 1.022726 1.022726
#> t9  0.000000 0.000000 0.000000 1.449474 1.460181 2.470506 1.022726 1.022726
#> t6  0.000000 0.000000 0.000000 1.022726 1.022726 1.022726 2.470506 1.345526
#> t1  0.000000 0.000000 0.000000 1.022726 1.022726 1.022726 1.345526 2.470506
#> t3  0.000000 0.000000 0.000000 1.022726 1.022726 1.022726 1.345526 1.413873
#> t2  0.000000 0.000000 0.000000 1.022726 1.022726 1.022726 1.343042 1.343042
#>           t3       t2
#> t4  0.000000 0.000000
#> t7  0.000000 0.000000
#> t5  0.000000 0.000000
#> t8  1.022726 1.022726
#> t10 1.022726 1.022726
#> t9  1.022726 1.022726
#> t6  1.345526 1.343042
#> t1  1.413873 1.343042
#> t3  2.470506 1.343042
#> t2  1.343042 2.470506
```
