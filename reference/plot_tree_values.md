# Plot a tree with branches coloured by a value

Colours every branch on a diverging scale (blue below `center`, red
above), e.g. by estimated slopes from
[`varying_slopes()`](https://thomaspueschel.github.io/phylogp/reference/varying_slopes.md)
or by ancestral states.

## Usage

``` r
plot_tree_values(
  tree,
  values,
  center = stats::median(values),
  lim = NULL,
  type = "phylogram",
  show.tip.label = NULL,
  edge.width = 2.5,
  main = NULL,
  legend = TRUE,
  key_labels = NULL,
  ...
)
```

## Arguments

- tree:

  A `phylo` tree.

- values:

  One value per edge, in the order of `tree$edge`.

- center:

  Value shown in the neutral colour.

- lim:

  Half-width of the colour range around `center`; values beyond it get
  the extreme colours. Default: the largest deviation.

- type:

  Tree layout, passed to
  [`ape::plot.phylo()`](https://rdrr.io/pkg/ape/man/plot.phylo.html).

- show.tip.label:

  Show tip labels? Default: only for small trees.

- edge.width:

  Branch width.

- main:

  Title.

- legend:

  Draw a colour key?

- key_labels:

  Labels for the low end, centre and high end of the key.

- ...:

  Further arguments to
  [`ape::plot.phylo()`](https://rdrr.io/pkg/ape/man/plot.phylo.html).

## Value

`NULL`, invisibly.
