# Distance-based kernels

A kernel `sigma^2 * phi(D / ell)` on a distance matrix: patristic
distances on the tree (the default) or any other distance between
species, such as great-circle distances between range centroids for a
phylogeny + space model (`k_bm() + k_dist(geo) + k_white()`). The
component is dense, so models containing it use dense linear algebra
(O(N^3)).

## Usage

``` r
k_dist(
  D = NULL,
  phi = c("exponential", "matern32", "matern52", "sqexp"),
  lower = 0.001,
  upper = 10
)
```

## Arguments

- D:

  `NULL` for patristic distances on the tree, or a symmetric distance
  matrix with row and column names matching the tip labels.

- phi:

  `"exponential"`, `"matern32"`, `"matern52"` or `"sqexp"`.

- lower, upper:

  Bounds for the length-scale `ell`, as multiples of the largest
  distance.

## Value

A kernel component.

## Details

Validity depends on the distance and on `phi`. On an ultrametric tree
every non-increasing `phi` gives a valid kernel on patristic distances.
On trees with fossil tips or without a clock only completely monotone
functions such as the exponential are guaranteed valid;
squared-exponential and Matern kernels can fail. Use
[`kernel_check()`](https://thomaspueschel.github.io/phylogp/reference/kernel_check.md)
before fitting.

## References

Freckleton R.P., Jetz W. 2009. Space versus phylogeny: Disentangling
phylogenetic and spatial signals in comparative data. Proceedings of the
Royal Society B. 276:21–30.

Guzmán Naranjo M., Becker L. 2022. Statistical bias control in typology.
Linguistic Typology. 26:605–670.

Hartmann F., Jäger G. 2023. Gaussian process models for geographic
controls in phylogenetic trees. Open Research Europe. 3:57.

Schoenberg I.J. 1938. Metric spaces and positive definite functions.
Transactions of the American Mathematical Society. 44:522–536.

Gneiting T. 2013. Strictly and non-strictly positive definite functions
on spheres. Bernoulli. 19:1327–1349.

Rasmussen C.E., Williams C.K.I. 2006. Gaussian processes for machine
learning. Cambridge, MA: MIT Press.

## Examples

``` r
k_dist(phi = "exponential")
#> phylogp kernel: dist(exponential) 
#>   dist1      dist(exponential)            parameters: sigma, ell
```
