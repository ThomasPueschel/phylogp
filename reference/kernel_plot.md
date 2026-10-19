# Plot the correlations implied by a kernel

Draws the tree next to the correlation matrix the kernel implies between
its tips, with tips in the same order, so clades appear as blocks.

## Usage

``` r
kernel_plot(
  kernel,
  tree,
  pars = NULL,
  labels = NULL,
  main = NULL,
  show_tree = TRUE
)
```

## Arguments

- kernel:

  A kernel, e.g. `k_bm() + k_white()`.

- tree:

  A `phylo` tree (or the output of
  [`tree_info()`](https://thomaspueschel.github.io/phylogp/reference/tree_info.md)).

- pars:

  Optional parameter values, as a list named by component id (see
  `print(kernel)`), e.g.
  `list(lambda1 = list(sigma = 1, lambda = 0.3))`. Unspecified values
  take the defaults of each component. For components with latent rates
  ([`k_relaxed()`](https://thomaspueschel.github.io/phylogp/reference/k_relaxed.md),
  [`k_rate_time()`](https://thomaspueschel.github.io/phylogp/reference/k_rate_time.md)),
  a rate path is drawn from the prior unless supplied (`z` or
  `lograte`).

- labels:

  Show tip labels? By default only for trees with at most 40 tips.

- main:

  Title.

- show_tree:

  Draw the tree next to the matrix? With `FALSE` only the matrix is
  drawn, into the current plotting region (useful for panels).

## Value

The correlation matrix, invisibly.

## Examples

``` r
tr <- ape::rcoal(24)
kernel_plot(k_bm() + k_white(), tr)

kernel_plot(k_ou(), tr, pars = list(ou1 = list(alpha = 5)))
```
