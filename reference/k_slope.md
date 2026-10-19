# Phylogenetically varying slopes

`k_slope(x) * k` multiplies a kernel `k` by `x_i x_j`, the kernel of a
regression slope on `x` whose deviations from the mean slope evolve
along the tree under `k`:

## Usage

``` r
k_slope(x, center = TRUE, latent = TRUE)
```

## Arguments

- x:

  Name of a column of the data passed to
  [`phylo_gp()`](https://thomaspueschel.github.io/phylogp/reference/phylo_gp.md),
  or a numeric vector named by tip label.

- center:

  Centre `x` at its mean over the fitted species?

- latent:

  Sample the slope deviations at the nodes (linear time) rather than use
  the dense marginal kernel.

## Value

An object to be multiplied by one kernel component.

## Details

\$\$y_i = \beta_0 + (\beta_1 + v_i) x_i + \dots, \quad v \sim GP(0,
k).\$\$

- `k_slope(x) * k_bm()` (or another branch-wise kernel): the slope
  evolves along the tree. The covariance is `sigma^2 (x x') o K`, valid
  whenever `K` is (Schur product theorem).

- `k_slope(x) * k_white()`: independent species-specific slope
  deviations, i.e. residual variance that grows with `x^2`. This is the
  natural control for a varying-slope model: a varying slope also makes
  the variance of `y` grow with distance from the centre of `x`, and
  comparing the two kernels separates evolving slopes from
  heteroscedasticity.

The mean slope is the regression coefficient of `x` in the model
formula, so `x` should normally appear there as well. Because the kernel
depends on where `x` is centred, `x` is centred at its mean by default;
report the centring with the results.

With `latent = TRUE` (default) the slope deviations at the nodes of the
tree are sampled as latent variables and the remaining Gaussian
components are integrated out by pruning, so the cost per gradient stays
linear in the number of species; this needs `method = "sample"`. With
`latent = FALSE` the dense marginal kernel is used (O(N^3)), which also
allows maximum likelihood.

## Examples

``` r
k_bm() + k_slope("body") * k_bm() + k_white()      # varying slopes
#> phylogp kernel: bm() + slope(body) * bm() + white() 
#>   bm1        bm()                         parameters: sigma
#>   slopebm1   slope(body) * bm()           parameters: sigma, (latent values)
#>   white1     white()                      parameters: sigma
k_bm() + k_slope("body") * k_white() + k_white()   # heteroscedastic control
#> phylogp kernel: bm() + slope(body) * white() + white() 
#>   bm1        bm()                         parameters: sigma
#>   slopewhite1 slope(body) * white()        parameters: sigma
#>   white1     white()                      parameters: sigma
```
