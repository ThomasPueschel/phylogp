# Tree bookkeeping shared by the kernels, the Stan data and the R evaluators.

#' Pre-compute the tree quantities used by phylogp
#'
#' Edges are put in ape "postorder" (children before parents), which is the
#' order the O(N) pruning likelihood and leave-one-out recursions need.
#'
#' @param tree An object of class `phylo` with branch lengths.
#' @return A list of class `phylogp_tree`.
#' @export
tree_info <- function(tree) {
  if (inherits(tree, "phylogp_tree")) return(tree)
  if (!inherits(tree, "phylo")) stop("`tree` must be a 'phylo' object.", call. = FALSE)
  if (is.null(tree$edge.length)) stop("`tree` must have branch lengths.", call. = FALSE)
  if (any(tree$edge.length < 0)) stop("`tree` has negative branch lengths.", call. = FALSE)
  regime <- mapped_edge(tree)                # regime map of a phytools simmap
  o <- ape::reorder.phylo(tree, "postorder", index.only = TRUE)
  tree <- ape::reorder.phylo(tree, "postorder")
  if (!is.null(regime)) {
    regime <- regime[o, , drop = FALSE]
    tree$mapped.edge <- regime
    if (!is.null(tree$maps)) tree$maps <- tree$maps[o]
  }
  N <- ape::Ntip(tree)
  E <- nrow(tree$edge)
  h <- ape::node.depth.edgelength(tree)
  anc <- tree$edge[, 1]
  des <- tree$edge[, 2]
  structure(list(
    tree = tree, N = N, E = E, n_nodes = N + tree$Nnode, root = N + 1L,
    anc = anc, des = des, len = tree$edge.length,
    hs = h[anc], he = h[des], tip_height = h[seq_len(N)],
    Tmax = max(h[seq_len(N)]),
    is_tip = as.numeric(des <= N),
    tip_edge = match(seq_len(N), des),
    first_child = as.integer(!duplicated(anc)),
    ultrametric = ape::is.ultrametric(tree, tol = 1e-6),
    labels = tree$tip.label,
    regime = if (is.null(regime)) NULL else unname_rows(regime),
    cache = new.env(parent = emptyenv())   # B, C and D, computed on first use
  ), class = "phylogp_tree")
}

unname_rows <- function(M) {
  rownames(M) <- NULL
  M
}

#' @export
print.phylogp_tree <- function(x, ...) {
  cat("phylogp tree:", x$N, "tips,", x$E, "edges, height", signif(x$Tmax, 4),
      if (x$ultrametric) "(ultrametric)" else "(not ultrametric)", "\n")
  invisible(x)
}

# N x E tip-by-edge incidence matrix: B[i, e] = 1 if edge e lies on the path
# from the root to tip i.
tree_B <- function(ti) cached(ti, "B", function() {
  B <- matrix(0, ti$N, ti$E)
  parent_edge <- integer(ti$n_nodes)
  parent_edge[ti$des] <- seq_len(ti$E)
  for (i in seq_len(ti$N)) {
    node <- i
    while (node != ti$root) {
      e <- parent_edge[node]
      B[i, e] <- 1
      node <- ti$anc[e]
    }
  }
  B
})

# Shared path length from the root to the MRCA of each pair of tips (t_ij)
tree_C <- function(ti) cached(ti, "C", function() {
  C <- ape::vcv.phylo(ti$tree)
  unname(C[ti$labels, ti$labels])
})

# Patristic distances between tips (d_ij)
tree_D <- function(ti) cached(ti, "D", function() {
  D <- ape::cophenetic.phylo(ti$tree)
  unname(D[ti$labels, ti$labels])
})

# Compute a tree matrix once per tree_info() object
cached <- function(ti, key, f) {
  env <- ti$cache
  if (is.null(env)) return(f())
  if (is.null(env[[key]])) env[[key]] <- f()
  env[[key]]
}

# Data block shared by every generated Stan program
stan_tree_data <- function(ti) {
  list(E = ti$E, n_nodes = ti$n_nodes, root = ti$root,
       anc = ti$anc, des = ti$des, first_child = ti$first_child,
       tip_edge = ti$tip_edge, len = ti$len, hs = ti$hs, he = ti$he,
       is_tip = ti$is_tip, Tmax = ti$Tmax)
}
