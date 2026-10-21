# User-defined branch-wise kernels

Defines a new tree kernel by the variance each edge accumulates, written
as a Stan expression. Any kernel of this form is a valid covariance on
every tree as long as the edge variances are non-negative, and it keeps
the O(N) pruning likelihood and leave-one-out densities.

## Usage

``` r
k_edge(w, pars = list(), defaults = NULL, name = "edge")
```

## Arguments

- w:

  Stan expression for the edge variances (before scaling by `sigma^2`).

- pars:

  Named list of parameter bounds, each `c(lower, upper)`. Parameters get
  uniform priors on these bounds.

- defaults:

  Named list of parameter values used by
  [`kernel_plot()`](https://thomaspueschel.github.io/phylogp/reference/kernel_plot.md)
  and
  [`kernel_check()`](https://thomaspueschel.github.io/phylogp/reference/kernel_check.md);
  defaults to the midpoints of the bounds.

- name:

  Short name used in printed output and Stan variable names.

## Value

A kernel component.

## Details

The expression must return a vector with one element per edge and may
use the edge vectors `len` (length), `hs` and `he` (heights of the start
and end of the edge above the root), `is_tip` (1 for tip edges), the
tree height `Tmax`, and the parameters named in `pars`. Use Stan's
elementwise operators `.*` and `./` between vectors. The same expression
is evaluated in R (after translating `.*`, `./`, `pow()` and `square()`)
for
[`kernel_matrix()`](https://thomaspueschel.github.io/phylogp/reference/kernel_matrix.md),
[`kernel_check()`](https://thomaspueschel.github.io/phylogp/reference/kernel_check.md)
and
[`kernel_plot()`](https://thomaspueschel.github.io/phylogp/reference/kernel_plot.md).

## References

Hansen T.F., Martins E.P. 1996. Translating between microevolutionary
process and macroevolutionary patterns: The correlation structure of
interspecific data. Evolution. 50:1404–1417.

Ho L. si T., Ané C. 2014. A linear-time algorithm for Gaussian and
non-Gaussian trait evolution models. Systematic Biology. 63:397–408.

## Examples

``` r
# Pagel's kappa with a rate trend through time
k <- k_edge("pow(len, kappa) .* exp(beta * he)",
            pars = list(kappa = c(0, 1), beta = c(-5, 5)))
k
#> phylogp kernel: edge("pow(len, kappa) .* exp(beta * he)") 
#>   edge1      edge("pow(len, kappa) .* exp(beta * he)") parameters: sigma, kappa, beta
```
