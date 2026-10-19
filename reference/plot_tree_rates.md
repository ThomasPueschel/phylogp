# Plot a tree with branches coloured by rate

Plot a tree with branches coloured by rate

## Usage

``` r
plot_tree_rates(
  tree,
  rates,
  lim = NULL,
  main = "Relative rate of evolution",
  ...
)
```

## Arguments

- tree:

  A `phylo` tree.

- rates:

  Positive rates, one per edge in the order of `tree$edge`.

- lim:

  Colour range on the log2 scale (rates from `2^-lim` to `2^lim`); by
  default symmetric around 1 and wide enough for all rates.

- main:

  Title.

- ...:

  Further arguments to
  [`plot_tree_values()`](https://thomaspueschel.github.io/phylogp/reference/plot_tree_values.md).

## Value

`NULL`, invisibly.
