# Pre-compute the tree quantities used by phylogp

Edges are put in ape "postorder" (children before parents), which is the
order the O(N) pruning likelihood and leave-one-out recursions need.

## Usage

``` r
tree_info(tree)
```

## Arguments

- tree:

  An object of class `phylo` with branch lengths.

## Value

A list of class `phylogp_tree`.
