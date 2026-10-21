# Plot regression effects

Plots the expected response, \\x^\top\beta\\, against one predictor at a
time, with a posterior credible band, over the data. Other predictors
are held at their mean (numeric) or reference level (factors). The line
is the population-level relationship, averaged over the phylogenetic and
other kernel effects, as in a mixed model. For maximum-likelihood fits
the line is drawn without a band.

## Usage

``` r
plot_effects(
  fit,
  effect = NULL,
  probs = c(0.05, 0.95),
  points = TRUE,
  n = 100,
  ndraws = 1000,
  plot = TRUE,
  col = "#4F6D8F",
  ...
)
```

## Arguments

- fit:

  A `phylogp_fit`.

- effect:

  Names of predictors (variables in the formula) to plot. By default
  every predictor, one panel each.

- probs:

  Lower and upper quantiles of the band.

- points:

  Draw the data? These are the species values, or the rows for fits with
  several rows per species.

- n:

  Number of grid values for a numeric predictor.

- ndraws:

  Number of posterior draws used.

- plot:

  Draw the plot? If `FALSE`, only the values are returned.

- col:

  Colour of the line and band.

- ...:

  Further graphical parameters passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly, a data frame with the predictor (`effect`), its value
(`value` for numeric predictors, `level` for factors), the posterior
median of the expected response (`estimate`) and the band (`lower`,
`upper`).

## Details

Named differently from `brms::conditional_effects()` so that both
packages can be loaded together.

## References

Bürkner P.-C. 2017. Brms: An R package for Bayesian multilevel models
using Stan. Journal of Statistical Software. 80:1–28.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- phylo_gp(brain ~ body, d, tree, k_bm() + k_white())
plot_effects(fit)
plot_effects(fit, "body", probs = c(0.025, 0.975))
} # }
```
