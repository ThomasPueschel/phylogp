# phylogp (development version)

* `predict()` now predicts species that were not fitted for every kernel,
  including those with latent values (`k_relaxed()`, `k_rate_time()`, latent
  `k_slope()` terms and `k_white(df = )`). The latent values a new species
  needs are drawn from their conditional prior given the fitted ones (a
  Brownian bridge along the fitted branch it joins, Brownian motion along its
  own branch, the prior for a Student-t multiplier).
* Several rows per species can now be combined with known measurement error
  (`se`, one value per row): the species means are weighted by
  `1 / (within_sigma^2 + se^2)`, and the likelihood stays exact and linear.
* Several rows per species can now be combined with a varying slope
  (`k_slope(x) * k_bm()`, latent): the slope acts on each row's own `x`, and
  leave-one-out integrates the left-out species' slope over all its rows in
  closed form.
* `data-raw/capabilities.R` checks every kernel against every feature; the
  result is in the new "What works with what" article.
* Bug fix: `loo()` for fits across a sample of trees (`phylo_gp_trees()`)
  now runs PSIS-LOO on each tree and averages each species' left-out density
  over trees, so trees are weighted equally, as documented. It used to run a
  single PSIS on the draws of all trees together, which weights each tree by
  the inverse of its predictive density for the left-out species (a harmonic
  mean across trees) and so underestimated elpd whenever trees differed in how
  well they predicted a species. The per-tree results are in
  `attr(, "per_tree")`.
* The paper is out as a preprint (EcoEvoRxiv, doi:10.32942/X2D99V);
  `citation("phylogp")` now gives it, with the software DOI.
* References: the help pages of the kernel components and tools, the
  vignettes and the README now cite the work each model builds on (Brownian
  motion, the mixed model, Pagel's transformations, early burst, OU, relaxed
  rates, spatial kernels, varying slopes, within-species variation, kernel
  validity, leave-one-out and prediction).
* `plot_effects()` plots the expected response along each predictor, with a
  posterior credible band, over the data (like `brms::conditional_effects()`,
  named differently so both packages can be loaded together), and returns the
  values for custom plots.
* New vignette, "A tour of models": from Brownian motion and PGLS through
  the classical family, measurement error, several individuals per species,
  rates through time, phylogeny plus geography and outlying species, to
  model choice by leave-one-out and leave-clade-out cross-validation,
  prediction and ancestral states.
* Several rows per species: with `species = ` given, `phylo_gp()` accepts
  repeated measurements (e.g. individuals) and adds a within-species residual
  standard deviation, `within_sigma`. The kernel then describes the species
  and `k_white()` is the non-phylogenetic species effect, as in
  `(1 | species)` of a mixed model. The likelihood is exact and keeps its
  linear cost, because it factorises into the species means and the
  deviations from them; predictors may vary within species. Leave-one-out,
  `loo_refit()` and `loo_group()` leave out whole species, `predict()`
  includes the within-species variation in `type = "response"`, and
  `ancestral()` works as before. Not yet available with varying slopes or
  known measurement error (`se`).
* `make_stancode()`, `make_standata()` and `stancode()` are renamed
  `kernel_stancode()` and `kernel_standata()`, so that phylogp and brms can be
  loaded together without masking each other's functions.
  `kernel_stancode()` takes a kernel or a fitted model (replacing
  `stancode(fit)`).

# phylogp 0.2.0

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
