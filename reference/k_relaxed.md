# Relaxed (autocorrelated) rates

The log rate of evolution itself evolves by Brownian motion along the
tree (as in autocorrelated relaxed clocks). The rate of edge e is
`sigma^2 * exp((z_start + z_end) / 2)`, with `z` a Brownian motion on
the tree with variance `tau^2` per unit time and `z = 0` at the root.
`tau = 0` gives constant-rate Brownian motion. The component keeps the
branch-wise form, so the likelihood stays linear in the number of tips.

## Usage

``` r
k_relaxed(tau_scale = 1)
```

## Arguments

- tau_scale:

  Scale of the half-normal prior on `tau`.

## Value

A kernel component.

## Details

The node values of `z` are latent variables, so models with this
component can only be fitted by sampling. In leave-one-out densities,
the log rate at the left-out tip is integrated out against its
conditional prior by Gauss-Hermite quadrature.

## Examples

``` r
k_relaxed() + k_white()
#> phylogp kernel: relaxed() + white() 
#>   relaxed1   relaxed()                    parameters: sigma, tau, (latent values)
#>   white1     white()                      parameters: sigma
```
