# Combine kernels

Adding kernels describes independent processes acting together; the
covariance of the sum is the sum of the covariances.

## Usage

``` r
# S3 method for class 'phylogp_kernel'
e1 + e2

# S3 method for class 'phylogp_component'
e1 + e2
```

## Arguments

- e1, e2:

  Kernels or kernel components.

## Value

A kernel of class `phylogp_kernel`.

## References

Rasmussen C.E., Williams C.K.I. 2006. Gaussian processes for machine
learning. Cambridge, MA: MIT Press.

Duvenaud D., Lloyd J., Grosse R., Tenenbaum J., Ghahramani Z. 2013.
Structure discovery in nonparametric regression through compositional
kernel search. Proceedings of the 30th international conference on
machine learning. 28:1166–1174.

## Examples

``` r
k <- k_bm() + k_white()    # Brownian motion plus residual variance
k
#> phylogp kernel: bm() + white() 
#>   bm1        bm()                         parameters: sigma
#>   white1     white()                      parameters: sigma
```
