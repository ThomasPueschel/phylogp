# phylogp (development version)

* `kernel_equivalent()` checks whether one kernel can reproduce the tip
  distribution of another: for each parameter setting of the first it
  searches over all parameters of the second for the smallest
  Kullback-Leibler divergence, comparing contrasts when the mean is unknown.
  It flags best-fitting parameters that sit at their bounds.
* `phylo_gp_trees()` fits a model across a sample of trees and pools the
  posterior (or collects maximum-likelihood estimates), reporting how much of
  the uncertainty comes from the phylogeny; `loo()` works on the pooled fit.
* `loo_group()` does leave-clade-out (or any leave-group-out)
  cross-validation: refits without each group and scores the group's joint
  predictive density, in linear time for branch-wise kernels. Returns a
  `kfold` object for `loo::loo_compare()`.
* Tree matrices (incidence, shared path lengths, distances) are cached per
  `tree_info()` object.

# phylogp 0.2.0

* `predict()` gives posterior predictions for fitted species and for species
  in the tree without data (kriging on the tree).
* `ancestral()` reconstructs ancestral states in linear time; for OU it uses
  the OU process itself rather than its tip kernel.
* `loo_refit()` replaces the PSIS-LOO contribution of influential species
  (high Pareto k) by an exact refit.
* Varying slopes: `k_slope(x) * k_bm()` (or another branch-wise kernel) for
  slopes that evolve along the tree, fitted with latent node values in linear
  time or with the dense marginal kernel; `k_slope(x) * k_white()` as the
  heteroscedastic control. `varying_slopes()` summarises the slopes.
* `k_white(df = )`: Student-t residuals for data with outlying species.
* Leave-one-out integrates each species' own latent values (slopes, rates,
  variance multipliers); `k_relaxed()` can now be combined with dense
  components.
* `plot_tree_values()` colours a tree by any value.
* New vignette: "Slopes, prediction and robust residuals".

# phylogp 0.1.0

* First version.
* Kernel grammar: classical kernels (`k_bm()`, `k_lambda()`, `k_kappa()`,
  `k_delta()`, `k_eb()`, `k_ou()`, `k_white()`), extensions (`k_relaxed()`,
  `k_rate_time()`, `k_dist()`) and user-defined branch-wise kernels
  (`k_edge()`), combined with `+`.
* `phylo_gp()` writes, compiles and fits the Stan program by maximum
  likelihood or HMC, with an O(N) pruning likelihood and O(N) leave-one-out
  densities for branch-wise kernels.
* Design tools: `kernel_check()`, `kernel_plot()`, `kernel_matrix()`,
  `kernel_simulate()`, `kernel_loglik()`, `make_stancode()`.
* Rate summaries: `edge_rates()`, `plot_edge_rates()`, `plot_tree_rates()`,
  `rate_through_time()`.
* Vignettes: "Getting started with phylogp" and "Designing kernels".
