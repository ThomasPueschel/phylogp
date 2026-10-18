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
