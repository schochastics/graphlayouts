# resolve edge weights into branch lengths, following the same NA/NULL/vector
# semantics as the stress layouts (see .layout_with_stress_dim)
.resolve_branch_lengths <- function(g, weights) {
  m <- igraph::ecount(g)
  if (is.null(weights)) {
    if ("weight" %in% igraph::edge_attr_names(g)) {
      return(as.numeric(igraph::E(g)$weight))
    }
    return(rep(1, m))
  }
  if (length(weights) == 1 && is.na(weights)) {
    return(rep(1, m))
  }
  if (length(weights) != m) {
    stop("weights must be a numeric vector of length ecount(g)", call. = FALSE)
  }
  as.numeric(weights)
}

# compute coordinates for a single tree (connected, acyclic subgraph)
.tree_unrooted_component <- function(
  sg,
  w_sub,
  mode,
  daylight_iter,
  iter,
  tol
) {
  n <- igraph::vcount(sg)
  if (n == 1) {
    return(matrix(c(0, 0), 1, 2))
  }
  if (n == 2) {
    return(matrix(c(0, 0, w_sub[1], 0), 2, 2, byrow = TRUE))
  }

  if (mode == "stress") {
    # stress majorization on patristic (tree path-length) distances
    D <- igraph::distances(sg, weights = w_sub)
    W <- 1 / D^2
    diag(W) <- 0
    xinit <- .init_layout(sg, D, mds = TRUE, n = n, dim = 2)
    return(stress_major(xinit, W, D, iter, tol))
  }

  # equal-angle / equal-daylight: build a rooted traversal
  root <- which.max(igraph::degree(sg))[1]
  res <- igraph::bfs(sg, root = root, parent = TRUE, order = TRUE)
  ord <- as.integer(res$order)
  father <- as.integer(res$parent)

  # subtree leaf counts (accumulate over reverse BFS order)
  is_parent <- logical(n)
  is_parent[father[!is.na(father)]] <- TRUE
  leaf_count <- integer(n)
  leaf_count[!is_parent] <- 1L
  for (v in rev(ord)) {
    f <- father[v]
    if (!is.na(f)) {
      leaf_count[f] <- leaf_count[f] + leaf_count[v]
    }
  }
  nleaves <- leaf_count[root]

  # branch length from each node to its parent (root: 0)
  branch_len <- numeric(n)
  non_root <- which(!is.na(father))
  if (length(non_root) > 0) {
    eids <- igraph::get_edge_ids(
      sg,
      as.vector(rbind(father[non_root], non_root))
    )
    branch_len[non_root] <- w_sub[eids]
  }

  parent0 <- father - 1L
  parent0[is.na(parent0)] <- -1L

  xy <- equal_angle_layout(
    ord - 1L,
    parent0,
    leaf_count,
    branch_len,
    nleaves
  )

  if (mode == "equaldaylight") {
    al <- igraph::as_adj_list(sg, mode = "all")
    adj_ptr <- integer(n + 1)
    for (i in seq_len(n)) {
      adj_ptr[i + 1] <- adj_ptr[i] + length(al[[i]])
    }
    adj_idx <- as.integer(unlist(al)) - 1L
    xy <- equal_daylight_layout(adj_ptr, adj_idx, ord - 1L, xy, daylight_iter)
  }

  xy
}

#' unrooted tree layout
#'
#' @name layout_tree_unrooted
#' @description arranges the nodes of a tree (or forest) in the plane without a
#' designated root, using algorithms developed for drawing phylogenies.
#' @param g igraph object. Each connected component must be a tree (acyclic).
#' @param mode which algorithm to use. One of `"equalangle"` (default),
#' `"equaldaylight"` or `"stress"`. See details.
#' @param weights possibly a numeric vector with edge weights, interpreted as
#' branch lengths. If this is NULL and the graph has a weight edge attribute,
#' then the attribute is used. If this is NA then unit branch lengths are used
#' (even if the graph has a weight attribute). By default, unit branch lengths
#' are used.
#' @param daylight_iter number of refinement sweeps for `mode = "equaldaylight"`.
#' @param iter number of iterations during stress optimization (only used for
#' `mode = "stress"`).
#' @param tol stopping criterion for stress optimization (only used for
#' `mode = "stress"`).
#' @param bbox width of layout. Only relevant to determine the placement of
#' disconnected components (forests).
#' @details Three algorithms are available:
#'
#' * `"equalangle"` (Felsenstein, 1989): each subtree is allotted an angular
#'   wedge proportional to its number of leaves. Linear time, but can look
#'   crowded on unbalanced trees.
#' * `"equaldaylight"` (Felsenstein): iteratively equalizes the angular gaps
#'   ("daylight") around each internal node, starting from the equal-angle
#'   layout, for a more balanced drawing.
#' * `"stress"`: stress majorization using the patristic (tree path-length)
#'   distances as target distances, reusing the package's stress engine.
#'
#' Branch lengths (edge weights) are respected when present; otherwise unit
#' lengths are used.
#'
#' The layout_igraph_* function should not be used directly. It is only used as
#' an argument for plotting with 'igraph'. 'ggraph' natively supports the layout.
#' @return matrix of xy coordinates
#' @references Felsenstein, J. (1989). PHYLIP - Phylogeny Inference Package
#' (Version 3.2). *Cladistics*, 5, 164-166.
#'
#' Bachmaier, C., Brandes, U., & Schlieper, B. (2005). Drawing phylogenetic
#' trees. *In International Symposium on Algorithms and Computation* (pp.
#' 1110-1121). Springer, Berlin, Heidelberg.
#' @examples
#' library(igraph)
#' g <- make_tree(20, 3, mode = "undirected")
#'
#' xy <- layout_as_tree_unrooted(g)
#' xy <- layout_as_tree_unrooted(g, mode = "equaldaylight")
#' @export
layout_as_tree_unrooted <- function(
  g,
  mode = c("equalangle", "equaldaylight", "stress"),
  weights = NA,
  daylight_iter = 15,
  iter = 500,
  tol = 0.0001,
  bbox = 30
) {
  ensure_igraph(g)
  mode <- match.arg(mode)

  oldseed <- get_seed()
  set.seed(42)
  on.exit(restore_seed(oldseed))

  w_full <- .resolve_branch_lengths(g, weights)

  comps <- igraph::components(g, "weak")

  # tag edges so weights can be mapped onto induced subgraphs
  igraph::edge_attr(g, "_edgename") <- seq_len(igraph::ecount(g))
  names(w_full) <- seq_len(igraph::ecount(g))

  lg <- list()
  node_order <- integer(0)

  for (i in seq_len(comps$no)) {
    idx <- which(comps$membership == i)
    sg <- igraph::induced_subgraph(g, idx)
    node_order <- c(node_order, idx)

    if (igraph::ecount(sg) != igraph::vcount(sg) - 1) {
      stop(
        "each connected component of g must be a tree (acyclic)",
        call. = FALSE
      )
    }

    edge_idx <- igraph::edge_attr(g, "_edgename") %in%
      igraph::edge_attr(sg, "_edgename")
    w_sub <- unname(w_full[edge_idx])

    lg[[i]] <- .tree_unrooted_component(
      sg,
      w_sub,
      mode,
      daylight_iter,
      iter,
      tol
    )
  }

  if (comps$no == 1) {
    return(lg[[1]])
  }

  lg <- lapply(lg, mv_to_null)
  p <- order(comps$csize)
  lg <- .component_mover(lg, p, bbox)
  x <- do.call("rbind", lg)
  x[order(node_order), , drop = FALSE]
}
