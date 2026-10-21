# Multiply a kernel by a covariate

See
[`k_slope()`](https://thomaspueschel.github.io/phylogp/reference/k_slope.md).

## Usage

``` r
# S3 method for class 'phylogp_scale'
e1 * e2

# S3 method for class 'phylogp_component'
e1 * e2

# S3 method for class 'phylogp_kernel'
e1 * e2
```

## Arguments

- e1, e2:

  A
  [`k_slope()`](https://thomaspueschel.github.io/phylogp/reference/k_slope.md)
  term and one kernel component.

## Value

A kernel component.

## References

Rasmussen C.E., Williams C.K.I. 2006. Gaussian processes for machine
learning. Cambridge, MA: MIT Press.

Duvenaud D., Lloyd J., Grosse R., Tenenbaum J., Ghahramani Z. 2013.
Structure discovery in nonparametric regression through compositional
kernel search. Proceedings of the 30th international conference on
machine learning. 28:1166–1174.

Hadfield J.D., Nakagawa S. 2010. General quantitative genetic methods
for comparative biology: Phylogenies, taxonomies and multi-trait models
for continuous and categorical characters. Journal of Evolutionary
Biology. 23:494–508.
