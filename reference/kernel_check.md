# Check that a kernel is a valid covariance on a tree

A kernel must be positive semi-definite on the tips of the tree. The
check first applies the theory (see the package vignette and Puschel, in
prep.):

## Usage

``` r
kernel_check(kernel, tree, n_grid = 25)
```

## Arguments

- kernel:

  A kernel, e.g. `k_bm() + k_white()`.

- tree:

  A `phylo` tree (or the output of
  [`tree_info()`](https://thomaspueschel.github.io/phylogp/reference/tree_info.md)).

- n_grid:

  Number of parameter values tried in numerical checks.

## Value

A data frame with one row per component (class `phylogp_check`), printed
as a short report.

## Details

- a branch-wise kernel (edge variances accumulated along the tree) is
  valid whenever its edge variances are non-negative;

- on an ultrametric tree, any non-increasing function of patristic
  distance is valid;

- on any tree, completely monotone functions of patristic distance, such
  as the exponential (Ornstein-Uhlenbeck) kernel, are valid;

- squared-exponential and Matern kernels on patristic distances are not
  guaranteed valid on trees with fossil tips or without a clock.

Where validity is not guaranteed, the check evaluates the kernel over a
grid of parameter values (edge variances for user-defined branch-wise
kernels, the smallest eigenvalue for dense kernels) and reports
failures. A sum of valid kernels is valid.

## References

Schoenberg I.J. 1938. Metric spaces and positive definite functions.
Transactions of the American Mathematical Society. 44:522–536.

Deza M.M., Laurent M. 1997. Geometry of cuts and metrics. Berlin:
Springer.

Dellacherie C., Martínez S., San Martín J. 2014. Inverse M-matrices and
ultrametric matrices. Cham: Springer.

Feragen A., Lauze F., Hauberg S. 2015. Geodesic exponential kernels:
When curvature and linearity conflict. Proceedings of the IEEE
conference on computer vision and pattern recognition (CVPR).:3032–3042.

## Examples

``` r
tr <- ape::rtree(30)                       # not ultrametric
kernel_check(k_dist(phi = "sqexp"), tr)    # may fail
#> Kernel check on a 30-tip non-ultrametric tree
#>   [INVALID] dist1      dist(sqexp) (dense)
#>          sqexp kernel on patristic distances of a non-ultrametric tree is not guaranteed valid; FAILS at 10 of 25 parameter values tried
#> The kernel is NOT guaranteed to be a valid covariance on this tree.
#> Use a branch-wise kernel or a completely monotone function of distance
#> (e.g. phi = "exponential"); adding jitter would change the model.
kernel_check(k_dist(phi = "exponential"), tr)
#> Kernel check on a 30-tip non-ultrametric tree
#>   [valid  ] dist1      dist(exponential) (dense)
#>          completely monotone function of a tree metric (valid on every tree)
#> The kernel is a valid covariance on this tree.
```
