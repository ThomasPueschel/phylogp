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

## Examples

``` r
k <- k_bm() + k_white()    # Brownian motion plus residual variance
k
#> phylogp kernel: bm() + white() 
#>   bm1        bm()                         parameters: sigma
#>   white1     white()                      parameters: sigma
```
