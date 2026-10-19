# Classical comparative-method kernels

Each function returns one kernel component. Components are combined with
`+`, e.g. `k_bm() + k_white()`. Every component has its own scale
parameter `sigma`; its contribution to the covariance is `sigma^2` times
the base kernel below (`t_ij`: shared path length from the root to the
most recent common ancestor of tips i and j; `d_ij`: patristic distance;
`T`: tree height).

## Usage

``` r
k_bm()

k_lambda(lower = 1e-07, upper = 1)

k_kappa(lower = 1e-06, upper = 1)

k_delta(lower = 1e-05, upper = 3)

k_eb(lower = function(ti, cmp) -log(1e+06)/ti$Tmax, upper = 0)

k_ou(
  root = c("random", "fixed"),
  lower = function(ti, cmp) 1e-07/ti$Tmax,
  upper = function(ti, cmp) 50/ti$Tmax
)

k_white(df = NULL)
```

## Arguments

- lower, upper:

  Bounds of the kernel parameter. They may be numbers or functions of
  the tree information (see
  [`tree_info()`](https://thomaspueschel.github.io/phylogp/reference/tree_info.md));
  the defaults match those of `phylolm`.

- root:

  For `k_ou()`, `"random"` (stationary) or `"fixed"` root.

- df:

  For `k_white()`: degrees of freedom of a Student-t residual (heavy
  tails, robust to outlying species). `NULL` (default) gives Gaussian
  residuals. The t residual is a scale mixture: each species gets a
  latent variance multiplier with an inverse-gamma(df/2, df/2) prior, so
  the model needs `method = "sample"`.

## Value

A kernel component.

## Details

- `k_bm()`: Brownian motion, `t_ij`.

- `k_lambda()`: Pagel's lambda. Internal branches are multiplied by
  `lambda` and tip branches lengthened to keep tip heights fixed.

- `k_kappa()`: Pagel's kappa, every branch length raised to `kappa`.

- `k_delta()`: Pagel's delta, node heights transformed, giving
  `T^(1 - delta) t_ij^delta`.

- `k_eb()`: early burst / ACDC, rate `exp(rate * t)`, giving
  `(exp(rate * t_ij) - 1) / rate`.

- `k_ou()`: Ornstein-Uhlenbeck with selection strength `alpha`; random
  (stationary) root `exp(-alpha d_ij) / (2 alpha)`, or fixed root, which
  multiplies this by `1 - exp(-2 alpha t_ij)`. On ultrametric trees it
  is fitted in branch-wise form (O(N)); on other trees through its dense
  closed form.

- `k_white()`: independent tip variance (residual, within-species
  variation); `k_bm() + k_white()` is Pagel's lambda on an ultrametric
  tree.

## Examples

``` r
k_bm() + k_white()
#> phylogp kernel: bm() + white() 
#>   bm1        bm()                         parameters: sigma
#>   white1     white()                      parameters: sigma
k_ou(root = "fixed")
#> phylogp kernel: ou(root = "fixed") 
#>   ou1        ou(root = "fixed")           parameters: sigma, alpha
```
