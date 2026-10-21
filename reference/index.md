# Package index

## Fitting models

Fit a phylogenetic Gaussian process and compare models.

- [`phylo_gp()`](https://thomaspueschel.github.io/phylogp/reference/phylo_gp.md)
  : Fit a phylogenetic Gaussian-process model
- [`phylo_gp_trees()`](https://thomaspueschel.github.io/phylogp/reference/phylo_gp_trees.md)
  : Fit a model across a sample of trees
- [`loo(`*`<phylogp_fit>`*`)`](https://thomaspueschel.github.io/phylogp/reference/loo.phylogp_fit.md)
  [`loo(`*`<phylogp_trees>`*`)`](https://thomaspueschel.github.io/phylogp/reference/loo.phylogp_fit.md)
  : Leave-one-out cross-validation for a phylogenetic GP fit
- [`loo_refit()`](https://thomaspueschel.github.io/phylogp/reference/loo_refit.md)
  : Exact leave-one-out for species with high Pareto k
- [`loo_group()`](https://thomaspueschel.github.io/phylogp/reference/loo_group.md)
  : Leave-clade-out cross-validation

## Prediction and ancestral states

Regression effects, species without data and internal nodes.

- [`plot_effects()`](https://thomaspueschel.github.io/phylogp/reference/plot_effects.md)
  : Plot regression effects
- [`predict(`*`<phylogp_fit>`*`)`](https://thomaspueschel.github.io/phylogp/reference/predict.phylogp_fit.md)
  : Predict trait values for species
- [`ancestral()`](https://thomaspueschel.github.io/phylogp/reference/ancestral.md)
  : Ancestral states

## Classical kernels

Brownian motion, Pagel’s transformations, early burst,
Ornstein-Uhlenbeck and residual variance. Combine them with `+`.

- [`k_bm()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md)
  [`k_lambda()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md)
  [`k_kappa()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md)
  [`k_delta()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md)
  [`k_eb()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md)
  [`k_ou()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md)
  [`k_white()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md)
  : Classical comparative-method kernels
- [`` `+`( ``*`<phylogp_kernel>`*`)`](https://thomaspueschel.github.io/phylogp/reference/kernel-sum.md)
  [`` `+`( ``*`<phylogp_component>`*`)`](https://thomaspueschel.github.io/phylogp/reference/kernel-sum.md)
  : Combine kernels

## Kernel extensions

Kernels for variable rates, distances and your own models.

- [`k_relaxed()`](https://thomaspueschel.github.io/phylogp/reference/k_relaxed.md)
  : Relaxed (autocorrelated) rates
- [`k_rate_time()`](https://thomaspueschel.github.io/phylogp/reference/k_rate_time.md)
  : Rate of evolution as a smooth function of time
- [`k_dist()`](https://thomaspueschel.github.io/phylogp/reference/k_dist.md)
  : Distance-based kernels
- [`k_edge()`](https://thomaspueschel.github.io/phylogp/reference/k_edge.md)
  : User-defined branch-wise kernels
- [`k_slope()`](https://thomaspueschel.github.io/phylogp/reference/k_slope.md)
  : Phylogenetically varying slopes
- [`` `*`( ``*`<phylogp_scale>`*`)`](https://thomaspueschel.github.io/phylogp/reference/kernel-product.md)
  [`` `*`( ``*`<phylogp_component>`*`)`](https://thomaspueschel.github.io/phylogp/reference/kernel-product.md)
  [`` `*`( ``*`<phylogp_kernel>`*`)`](https://thomaspueschel.github.io/phylogp/reference/kernel-product.md)
  : Multiply a kernel by a covariate

## Designing and checking kernels

Look at a kernel, simulate from it and check that it is valid on a tree.

- [`kernel_check()`](https://thomaspueschel.github.io/phylogp/reference/kernel_check.md)
  : Check that a kernel is a valid covariance on a tree
- [`kernel_plot()`](https://thomaspueschel.github.io/phylogp/reference/kernel_plot.md)
  : Plot the correlations implied by a kernel
- [`kernel_matrix()`](https://thomaspueschel.github.io/phylogp/reference/kernel_matrix.md)
  : Covariance matrix implied by a kernel
- [`kernel_simulate()`](https://thomaspueschel.github.io/phylogp/reference/kernel_simulate.md)
  : Simulate traits from a kernel
- [`kernel_loglik()`](https://thomaspueschel.github.io/phylogp/reference/kernel_loglik.md)
  : Log-likelihood of a kernel in R
- [`kernel_equivalent()`](https://thomaspueschel.github.io/phylogp/reference/kernel_equivalent.md)
  : Do two kernels imply the same distribution for the tips?

## Rates, slopes and tree plots

Summaries of estimated rates and slopes, and trees coloured by them.

- [`edge_rates()`](https://thomaspueschel.github.io/phylogp/reference/edge_rates.md)
  : Branch rates from a relaxed-rate fit
- [`plot_edge_rates()`](https://thomaspueschel.github.io/phylogp/reference/plot_edge_rates.md)
  : Plot branch rates on the tree
- [`plot_tree_rates()`](https://thomaspueschel.github.io/phylogp/reference/plot_tree_rates.md)
  : Plot a tree with branches coloured by rate
- [`rate_through_time()`](https://thomaspueschel.github.io/phylogp/reference/rate_through_time.md)
  : Rate of evolution through time
- [`varying_slopes()`](https://thomaspueschel.github.io/phylogp/reference/varying_slopes.md)
  : Species- and node-level slopes from a varying-slope fit
- [`plot_tree_values()`](https://thomaspueschel.github.io/phylogp/reference/plot_tree_values.md)
  : Plot a tree with branches coloured by a value

## Stan code and data

- [`kernel_stancode()`](https://thomaspueschel.github.io/phylogp/reference/kernel_stancode.md)
  : Stan code for a kernel or a fitted model
- [`kernel_standata()`](https://thomaspueschel.github.io/phylogp/reference/kernel_standata.md)
  : Data for the Stan program of a kernel
- [`tree_info()`](https://thomaspueschel.github.io/phylogp/reference/tree_info.md)
  : Pre-compute the tree quantities used by phylogp
