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

## References

Thorne J.L., Kishino H. 2002. Divergence time and evolutionary rate
estimation with multilocus data. Systematic Biology. 51:689–702.

O’Meara B.C., Ané C., Sanderson M.J., Wainwright P.C. 2006. Testing for
different rates of continuous trait evolution using likelihood.
Evolution. 60:922–933.

Eastman J.M., Alfaro M.E., Joyce P., Hipp A.L., Harmon L.J. 2011. A
novel comparative method for identifying shifts in the rate of character
evolution on trees. Evolution. 65:3578–3589.

Venditti C., Meade A., Pagel M. 2011. Multiple routes to mammalian
diversity. Nature. 479:393–396.

Revell L.J. 2021. A variable-rate quantitative trait evolution model
using penalized-likelihood. PeerJ. 9:e11997.

Gill M.S., Ho L.S.T., Baele G., Lemey P., Suchard M.A. 2017. A relaxed
directional random walk model for phylogenetic trait evolution.
Systematic Biology. 66:299–319.

Fisher A.A., Ji X., Zhang Z., Lemey P., Suchard M.A. 2021. Relaxed
random walks at scale. Systematic Biology. 70:258–267.

## Examples

``` r
k_relaxed() + k_white()
#> phylogp kernel: relaxed() + white() 
#>   relaxed1   relaxed()                    parameters: sigma, tau, (latent values)
#>   white1     white()                      parameters: sigma
```
