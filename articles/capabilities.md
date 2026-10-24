# What works with what

Each kernel was fitted to a small simulated data set (24 species, one in
the tree but not fitted) and every feature was tried; ✓ means it ran, –
means it is not defined for that kernel (see the notes). Kernels without
a
[`k_white()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md)
term behave like the same kernel with one. The table is made by
`data-raw/capabilities.R` in the package sources, so it can be rerun
after any change.

| Kernel | Maximum likelihood | Bayesian (HMC) | loo() | loo_refit() | loo_group() | predict(), fitted species | predict(), new species | ancestral() | se | Several rows per species | … loo() | … predict(), new species | … with se | phylo_gp_trees() + loo() |
|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|
| `k_bm() + k_white()` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| [`k_lambda()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `k_kappa() + k_white()` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `k_delta() + k_white()` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `k_eb() + k_white()` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `k_ou() + k_white()` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `k_bm() + k_white(df = 4)` | – ^(a) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `k_relaxed() + k_white()` | – ^(a) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `k_rate_time(5) + k_white()` | – ^(a) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `k_bm() + k_dist() + k_white()` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – ^(b) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `k_edge() + k_white()` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `k_bm() + k_slope(x) * k_bm() + k_white()` | – ^(a) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – ^(c) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `k_bm() + k_slope(x, latent = FALSE) * k_bm() + k_white()` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – ^(b) | ✓ | – ^(d) |  |  | – ^(d) | ✓ |
| `k_bm() + k_slope(x) * k_white()` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – ^(d) |  |  | – ^(d) | ✓ |

Columns starting with … refer to several rows per species.

^(a) Latent values (rates, slopes, t multipliers) are sampled, so these
kernels need `method = "sample"`.  
^(b) Ancestral states need a process along the tree; dense components
([`k_dist()`](https://thomaspueschel.github.io/phylogp/reference/k_dist.md),
`k_slope(latent = FALSE)`) have none.  
^(c) Ancestral states are not defined for varying slopes (the slope at a
node has no trait value).  
^(d) With several rows per species, use one latent `k_slope(x) * k` term
with a branch-wise `k`.
