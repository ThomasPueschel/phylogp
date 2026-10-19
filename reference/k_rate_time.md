# Rate of evolution as a smooth function of time

The rate `sigma^2(t)` is piecewise constant over `n_bins` equal time
bins from the root to the tips, with a Gaussian random-walk prior on its
logarithm. Each edge accumulates the integral of the rate over its time
span, so the component keeps the branch-wise form. Early burst is a
special case (a linear trend in log rate).

## Usage

``` r
k_rate_time(n_bins = 20, tau_scale = 0.5)
```

## Arguments

- n_bins:

  Number of time bins.

- tau_scale:

  Scale of the half-normal prior on the random-walk standard deviation
  `tau`.

## Value

A kernel component.

## Details

The bin values are latent variables, so models with this component can
only be fitted by sampling.

## Examples

``` r
k_rate_time(n_bins = 20) + k_white()
#> phylogp kernel: rate_time(n_bins = 20) + white() 
#>   ratetime1  rate_time(n_bins = 20)       parameters: sigma, tau, (latent values)
#>   white1     white()                      parameters: sigma
```
