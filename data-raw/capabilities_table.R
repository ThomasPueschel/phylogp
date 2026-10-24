# Writes the "What works with what" article from data-raw/capabilities.csv
# (made by data-raw/capabilities.R). Run from the package root.

x <- read.csv("data-raw/capabilities.csv", check.names = FALSE)
cols <- c(ml = "Maximum likelihood", hmc = "Bayesian (HMC)", loo = "loo()",
          loo_refit = "loo_refit()", loo_group = "loo_group()",
          predict_fitted = "predict(), fitted species", predict_new = "predict(), new species",
          ancestral = "ancestral()", se = "se", repeated = "Several rows per species",
          repeated_loo = "... loo()", repeated_predict_new = "... predict(), new species",
          repeated_se = "... with se", trees = "phylo_gp_trees() + loo()")
why <- c("need method = \"sample\"" = "a",
         "Ancestral states need a kernel that is a process along the tree" = "b",
         "Ancestral states are not available for models with varying slopes" = "c",
         "varying slopes are available as one" = "d")
notes <- c(a = "Latent values (rates, slopes, t multipliers) are sampled, so these kernels need `method = \"sample\"`.",
           b = "Ancestral states need a process along the tree; dense components (`k_dist()`, `k_slope(latent = FALSE)`) have none.",
           c = "Ancestral states are not defined for varying slopes (the slope at a node has no trait value).",
           d = "With several rows per species, use one latent `k_slope(x) * k` term with a branch-wise `k`.")
cell <- function(v) {
  if (is.na(v)) return("")
  if (v == "ok") return("✓")
  k <- names(why)[vapply(names(why), grepl, NA, x = v, fixed = TRUE)]
  if (length(k)) sprintf("– <sup>%s</sup>", why[k[1]]) else "✗"
}
tab <- vapply(names(cols), function(cl) vapply(x[[cl]], cell, ""), character(nrow(x)))
tab <- matrix(tab, nrow(x))
md <- c(paste0("| Kernel | ", paste(cols, collapse = " | "), " |"),
        paste0("|", paste(rep("---", length(cols) + 1), collapse = "|"), "|"),
        paste0("| `", x$kernel, "` | ", apply(tab, 1, paste, collapse = " | "), " |"))
if (any(tab == "✗")) warning("Unexplained failures: check data-raw/capabilities.csv")

writeLines(c(
  "---",
  "title: \"What works with what\"",
  "---",
  "",
  "Each kernel was fitted to a small simulated data set (24 species, one in the",
  "tree but not fitted) and every feature was tried; ✓ means it ran, – means it",
  "is not defined for that kernel (see the notes). Kernels without a `k_white()`",
  "term behave like the same kernel with one. The table is made by",
  "`data-raw/capabilities.R` in the package sources, so it can be rerun after any",
  "change.",
  "",
  md,
  "",
  "Columns starting with ... refer to several rows per species.",
  "",
  paste0("<sup>", names(notes), "</sup> ", notes, collapse = "  \n")),
  "vignettes/articles/capabilities.Rmd")
