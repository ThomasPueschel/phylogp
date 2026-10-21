# Changelog

## phylogp (development version)

- Several rows per species: with `species =` given,
  [`phylo_gp()`](https://thomaspueschel.github.io/phylogp/reference/phylo_gp.md)
  accepts repeated measurements (e.g. individuals) and adds a
  within-species residual standard deviation, `within_sigma`. The kernel
  then describes the species and
  [`k_white()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md)
  is the non-phylogenetic species effect, as in `(1 | species)` of a
  mixed model. The likelihood is exact and keeps its linear cost,
  because it factorises into the species means and the deviations from
  them; predictors may vary within species. Leave-one-out,
  [`loo_refit()`](https://thomaspueschel.github.io/phylogp/reference/loo_refit.md)
  and
  [`loo_group()`](https://thomaspueschel.github.io/phylogp/reference/loo_group.md)
  leave out whole species,
  [`predict()`](https://rdrr.io/r/stats/predict.html) includes the
  within-species variation in `type = "response"`, and
  [`ancestral()`](https://thomaspueschel.github.io/phylogp/reference/ancestral.md)
  works as before. Not yet available with varying slopes or known
  measurement error (`se`).
- `make_stancode()`, `make_standata()` and `stancode()` are renamed
  [`kernel_stancode()`](https://thomaspueschel.github.io/phylogp/reference/kernel_stancode.md)
  and
  [`kernel_standata()`](https://thomaspueschel.github.io/phylogp/reference/kernel_standata.md),
  so that phylogp and brms can be loaded together without masking each
  other’s functions.
  [`kernel_stancode()`](https://thomaspueschel.github.io/phylogp/reference/kernel_stancode.md)
  takes a kernel or a fitted model (replacing `stancode(fit)`).

## phylogp 0.2.0

- [`kernel_equivalent()`](https://thomaspueschel.github.io/phylogp/reference/kernel_equivalent.md)
  checks whether one kernel can reproduce the tip distribution of
  another: for each parameter setting of the first it searches over all
  parameters of the second for the smallest Kullback-Leibler divergence,
  comparing contrasts when the mean is unknown. It flags best-fitting
  parameters that sit at their bounds.

- [`phylo_gp_trees()`](https://thomaspueschel.github.io/phylogp/reference/phylo_gp_trees.md)
  fits a model across a sample of trees and pools the posterior (or
  collects maximum-likelihood estimates), reporting how much of the
  uncertainty comes from the phylogeny;
  [`loo()`](https://mc-stan.org/loo/reference/loo.html) works on the
  pooled fit.

- [`loo_group()`](https://thomaspueschel.github.io/phylogp/reference/loo_group.md)
  does leave-clade-out (or any leave-group-out) cross-validation: refits
  without each group and scores the group’s joint predictive density, in
  linear time for branch-wise kernels. Returns a `kfold` object for
  [`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html).

- Tree matrices (incidence, shared path lengths, distances) are cached
  per
  [`tree_info()`](https://thomaspueschel.github.io/phylogp/reference/tree_info.md)
  object.

- [`predict()`](https://rdrr.io/r/stats/predict.html) gives posterior
  predictions for fitted species and for species in the tree without
  data (kriging on the tree).

- [`ancestral()`](https://thomaspueschel.github.io/phylogp/reference/ancestral.md)
  reconstructs ancestral states in linear time; for OU it uses the OU
  process itself rather than its tip kernel.

- [`loo_refit()`](https://thomaspueschel.github.io/phylogp/reference/loo_refit.md)
  replaces the PSIS-LOO contribution of influential species (high
  Pareto k) by an exact refit.

- Varying slopes: `k_slope(x) * k_bm()` (or another branch-wise kernel)
  for slopes that evolve along the tree, fitted with latent node values
  in linear time or with the dense marginal kernel;
  `k_slope(x) * k_white()` as the heteroscedastic control.
  [`varying_slopes()`](https://thomaspueschel.github.io/phylogp/reference/varying_slopes.md)
  summarises the slopes.

- `k_white(df = )`: Student-t residuals for data with outlying species.

- Leave-one-out integrates each species’ own latent values (slopes,
  rates, variance multipliers);
  [`k_relaxed()`](https://thomaspueschel.github.io/phylogp/reference/k_relaxed.md)
  can now be combined with dense components.

- [`plot_tree_values()`](https://thomaspueschel.github.io/phylogp/reference/plot_tree_values.md)
  colours a tree by any value.

- New vignette: “Slopes, prediction and robust residuals”.

## phylogp 0.1.0

- First version.
- Kernel grammar: classical kernels
  ([`k_bm()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md),
  [`k_lambda()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md),
  [`k_kappa()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md),
  [`k_delta()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md),
  [`k_eb()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md),
  [`k_ou()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md),
  [`k_white()`](https://thomaspueschel.github.io/phylogp/reference/classical-kernels.md)),
  extensions
  ([`k_relaxed()`](https://thomaspueschel.github.io/phylogp/reference/k_relaxed.md),
  [`k_rate_time()`](https://thomaspueschel.github.io/phylogp/reference/k_rate_time.md),
  [`k_dist()`](https://thomaspueschel.github.io/phylogp/reference/k_dist.md))
  and user-defined branch-wise kernels
  ([`k_edge()`](https://thomaspueschel.github.io/phylogp/reference/k_edge.md)),
  combined with `+`.
- [`phylo_gp()`](https://thomaspueschel.github.io/phylogp/reference/phylo_gp.md)
  writes, compiles and fits the Stan program by maximum likelihood or
  HMC, with an O(N) pruning likelihood and O(N) leave-one-out densities
  for branch-wise kernels.
- Design tools:
  [`kernel_check()`](https://thomaspueschel.github.io/phylogp/reference/kernel_check.md),
  [`kernel_plot()`](https://thomaspueschel.github.io/phylogp/reference/kernel_plot.md),
  [`kernel_matrix()`](https://thomaspueschel.github.io/phylogp/reference/kernel_matrix.md),
  [`kernel_simulate()`](https://thomaspueschel.github.io/phylogp/reference/kernel_simulate.md),
  [`kernel_loglik()`](https://thomaspueschel.github.io/phylogp/reference/kernel_loglik.md),
  `make_stancode()`.
- Rate summaries:
  [`edge_rates()`](https://thomaspueschel.github.io/phylogp/reference/edge_rates.md),
  [`plot_edge_rates()`](https://thomaspueschel.github.io/phylogp/reference/plot_edge_rates.md),
  [`plot_tree_rates()`](https://thomaspueschel.github.io/phylogp/reference/plot_tree_rates.md),
  [`rate_through_time()`](https://thomaspueschel.github.io/phylogp/reference/rate_through_time.md).
- Vignettes: “Getting started with phylogp” and “Designing kernels”.
