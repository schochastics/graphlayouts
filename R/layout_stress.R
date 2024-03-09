.init_layout <- function(g, D, mds, n, dim) {
    if (!mds) {
        return(matrix(stats::runif(n * dim, 0, 1), n, dim))
    } else {
        rmat <- matrix(stats::runif(n * dim, -0.1, 0.1), n, dim)
        if (n <= 100) {
            return(igraph::layout_with_mds(g, dim = dim) + rmat)
        } else {
            return(layout_with_pmds(g, D = D[, sample(1:n, 100)], dim = dim) + rmat)
        }
    }
}

.component_mover <- function(lg, p, bbox) {
    curx <- 0
    cury <- 0
    maxy <- 0
    for (comp in p) {
        if (curx + max(lg[[comp]][, 1]) > bbox) {
            curx <- 0
            cury <- maxy + 1
        }
        lg[[comp]][, 1] <- lg[[comp]][, 1] + curx
        lg[[comp]][, 2] <- lg[[comp]][, 2] + cury
        curx <- max(lg[[comp]][, 1]) + 1
        maxy <- max(c(maxy, max(lg[[comp]][, 2])))
    }
    return(lg)
}

.component_layouter <- function(g, weights, comps, dim, mds, bbox, iter, tol, FUN, ...) {
    # check which ... are arguments of FUN
    FUN <- match.fun(FUN)
    params_in <- list(...)
    FUN_formals <- formals(FUN)
    idx <- names(params_in) %in% names(FUN_formals)
    params <- params_in[idx]
    if ("dim" %in% names(FUN_formals)) {
        params <- c(params, list(dim = params_in[["fixdim"]]))
    }

    lg <- list()
    node_order <- c()
    if (!is.null(weights) && any(!is.na(weights))) {
        igraph::edge_attr(g, "_edgename") <- 1:igraph::ecount(g)
        names(weights) <- 1:igraph::ecount(g)
    }

    for (i in 1:comps$no) {
        idx <- which(comps$membership == i)
        sg <- igraph::induced_subgraph(g, idx)
        edge_idx <- igraph::edge_attr(g, "_edgename") %in% igraph::edge_attr(sg, "_edgename")
        n <- igraph::vcount(sg)
        node_order <- c(node_order, idx)

        if (n == 1) {
            lg[[i]] <- matrix(rep(0, dim), 1, dim, byrow = TRUE)
        } else if (n == 2) {
            lg[[i]] <- matrix(c(0, rep(0, dim - 1), 1, rep(0, dim - 1)), 2, dim, byrow = TRUE)
            next()
        } else {
            if (!is.null(weights) && any(!is.na(weights))) {
                D <- igraph::distances(sg, weights = weights[edge_idx])
            } else {
                D <- igraph::distances(sg, weights = weights)
            }
            W <- 1 / D^2
            diag(W) <- 0

            xinit <- .init_layout(sg, D, mds, n, dim)
            if ("dim" %in% names(params)) {
                xinit[, params[["dim"]]] <- params_in[["coord"]][idx]
            }
            params_FUN <- c(params, list(y = xinit, W = W, D = D, iter = iter, tol = tol))
            lg[[i]] <- do.call(FUN, params_FUN) # FUN(xinit, W, D, iter, tol)
        }
    }
    if (!"dim" %in% names(params)) {
        lg <- lapply(lg, mv_to_null)
        p <- order(comps$csize)
        lg <- .component_mover(lg, p, bbox)
    }
    x <- do.call("rbind", lg)
    x[order(node_order), , drop = FALSE]
}

.layout_with_stress_dim <- function(g, weights = NA, iter = 500, tol = 0.0001, mds = TRUE, bbox = 30, dim = 2) {
    ensure_igraph(g)
    if (!dim %in% c(2, 3)) {
        stop("dim must be either 2 or 3")
    }

    oldseed <- get_seed()
    set.seed(42) # stress is deterministic and produces the same result up to translation. This keeps the layout fixed
    on.exit(restore_seed(oldseed))

    comps <- igraph::components(g, "weak")
    if (comps$no == 1) {
        n <- igraph::vcount(g)

        if (n == 1) {
            return(matrix(rep(0, dim), 1, dim, byrow = TRUE))
        } else if (n == 2) {
            return(matrix(c(0, rep(0, dim - 1), 1, rep(0, dim - 2)), 2, dim, byrow = TRUE))
        } else {
            if (!is.null(weights) && any(!is.na(weights))) {
                D <- igraph::distances(g, weights = weights)
            } else {
                D <- igraph::distances(g)
            }
            W <- 1 / D^2
            diag(W) <- 0

            xinit <- .init_layout(g, D, mds, n, dim)

            if (dim == 2) {
                return(stress_major(xinit, W, D, iter, tol))
            } else {
                return(stress_major3D(xinit, W, D, iter, tol))
            }
        }
    } else {
        layouter <- ifelse(dim == 2, stress_major, stress_major3D)
        return(.component_layouter(
            g = g, weights = weights, comps = comps, dim = dim, mds = mds,
            bbox = bbox, iter = iter, tol = tol, FUN = layouter
        ))
    }
}


#' stress majorization layout
#'
#' @name layout_stress
#' @rdname layout_stress
#' @description force-directed graph layout based on stress majorization.
#' Similar to Kamada-Kawai, but generally faster and with better results.
#' @param g igraph object
#' @param weights possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @param mds should an MDS layout be used as initial layout (default: TRUE)
#' @param bbox width of layout. Only relevant to determine the placement of disconnected graphs
#' @details Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight).
#'
#' The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' 'ggraph' natively supports the layout.
#' @seealso [layout_stress3D]
#' @return matrix of xy coordinates
#' @references Gansner, E. R., Koren, Y., & North, S. (2004). Graph drawing by stress majorization. *In International Symposium on Graph Drawing* (pp. 239-250). Springer, Berlin, Heidelberg.
#' @examples
#' library(igraph)
#' library(ggraph)
#' set.seed(665)
#'
#' g <- sample_pa(100, 1, 1, directed = FALSE)
#'
#' # calculate layout manually
#' xy <- layout_with_stress(g)
#'
#' # use it with ggraph
#' \dontrun{
#' ggraph(g, layout = "stress") +
#'     geom_edge_link0(edge_width = 0.2, colour = "grey") +
#'     geom_node_point(col = "black", size = 0.3) +
#'     theme_graph()
#' }
#' @export
layout_with_stress <- function(g, weights = NA, iter = 500, tol = 0.0001, mds = TRUE, bbox = 30) {
    .layout_with_stress_dim(g, weights, iter, tol, mds, bbox, dim = 2)
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' stress majorization layout in 3D
#'
#' @name layout_stress3D
#' @rdname layout_stress3D
#' @description force-directed graph layout based on stress majorization in 3D.
#' @param g igraph object
#' @param weights possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @param mds should an MDS layout be used as initial layout (default: TRUE)
#' @param bbox width of layout. Only relevant to determine the placement of disconnected graphs
#' @details Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight).
#'
#' @return matrix of xyz coordinates
#' @seealso [layout_stress]
#' @references Gansner, E. R., Koren, Y., & North, S. (2004). Graph drawing by stress majorization. *In International Symposium on Graph Drawing* (pp. 239-250). Springer, Berlin, Heidelberg.
#' @export
layout_with_stress3D <- function(g, weights = NA, iter = 500, tol = 0.0001, mds = TRUE, bbox = 30) {
    .layout_with_stress_dim(g, weights, iter, tol, mds, bbox, dim = 3)
}


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' radial focus layout
#'
#' @description arrange nodes in concentric circles around a focal node according to their distance from the focus.
#'
#' @name layout_focus
#' @param g igraph object
#' @param v id of focal node to be placed in the center
#' @param weights possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @details Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight).
#' @seealso [layout_focus_group]
#' The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' 'ggraph' natively supports the layout.
#' @return a list containing xy coordinates and the distances to the focal node
#' @references Brandes, U., & Pich, C. (2011). More flexible radial layout. *Journal of Graph Algorithms and Applications*, 15(1), 157-173.
#' @examples
#' library(igraph)
#' g <- sample_gnp(10, 0.4)
#' coords <- layout_with_focus(g, v = 1)
#' coords
#' @export

layout_with_focus <- function(g, v, weights = NA, iter = 500, tol = 0.0001) {
    ensure_igraph(g)
    ensure_connected(g)
    if (missing(v)) {
        stop("v missing without a default")
    }
    oldseed <- get_seed()
    set.seed(42) # stress is deterministic and produces same result up to translation. This keeps the layout fixed
    on.exit(restore_seed(oldseed))

    n <- igraph::vcount(g)
    D <- igraph::distances(g, weights = weights)
    W <- 1 / D^2
    diag(W) <- 0

    Z <- matrix(0, n, n)
    Z[v, ] <- Z[, v] <- 1
    Z <- W * Z


    rmat <- matrix(stats::runif(n * 2, -0.1, 0.1), n, 2)
    xinit <- igraph::layout_with_mds(g) + rmat

    tseq <- seq(0, 1, 0.1)
    x <- stress_focus(xinit, W, D, Z, tseq, iter, tol)

    offset <- x[v, ]
    x <- t(apply(x, 1, function(x) x - offset))
    return(list(xy = x, distance = D[, v]))
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' radial centrality layout
#'
#' @description arranges nodes in concentric circles according to a centrality index.
#'
#' @name layout_centrality
#' @param g igraph object
#' @param cent centrality scores
#' @param scale logical. should centrality scores be scaled to \eqn{[0,100]}? (Default: TRUE)
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @param tseq numeric vector. increasing sequence of coefficients to combine regular stress and constraint stress. See details.
#' @details The function optimizes a convex combination of regular stress and a constrained stress function which forces
#' nodes to be arranged on concentric circles. The vector `tseq` is the sequence of parameters used for the convex combination.
#' In iteration i of the algorithm \eqn{tseq[i]} is used to combine regular and constraint stress as \eqn{(1-tseq[i])*stress_{regular}+tseq[i]*stress_{constraint}}. The sequence must be increasing, start at zero and end at one. The default setting should be a good choice for most graphs.
#'
#' The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' 'ggraph' natively supports the layout.
#' @seealso [layout_centrality_group]
#' @return matrix of xy coordinates
#' @references Brandes, U., & Pich, C. (2011). More flexible radial layout. Journal of Graph Algorithms and Applications, 15(1), 157-173.
#' @examples
#' library(igraph)
#' library(ggraph)
#'
#' g <- sample_gnp(10, 0.4)
#' \dontrun{
#' ggraph(g, layout = "centrality", centrality = closeness(g)) +
#'     draw_circle(use = "cent") +
#'     geom_edge_link0() +
#'     geom_node_point(shape = 21, fill = "grey25", size = 5) +
#'     theme_graph() +
#'     coord_fixed()
#' }

#' @export
#'
layout_with_centrality <- function(g, cent, scale = TRUE, iter = 500, tol = 0.0001, tseq = seq(0, 1, 0.2)) {
    ensure_igraph(g)
    ensure_connected(g)
    if (missing(cent)) {
        stop("cent missing without a default")
    }
    oldseed <- get_seed()
    set.seed(42) # stress is deterministic and produces same result up to translation. This keeps the layout fixed
    on.exit(restore_seed(oldseed))

    n <- igraph::vcount(g)
    if (scale) {
        cent <- scale_to_100(cent)
    }
    r <- unname(igraph::diameter(g) / 2 * (1 - ((cent - min(cent)) / (max(cent) - min(cent) + 1))))

    D <- igraph::distances(g, weights = NA)
    W <- 1 / D^2
    diag(W) <- 0

    rmat <- matrix(stats::runif(n * 2, -0.1, 0.1), n, 2)
    xinit <- igraph::layout_with_mds(g) + rmat

    x <- stress_major(xinit, W, D, iter, tol)
    x <- stress_radii(x, W, D, r, tseq)

    # move highest cent to 0,0
    idx <- which.max(cent)[1]
    offset <- x[idx, ]

    x <- t(apply(x, 1, function(x) x - offset))
    if (scale) {
        radii_new <- round(100 - cent, 8)
        angles <- apply(x, 1, function(y) atan2(y[2], y[1]))
        return(cbind(radii_new * cos(angles), radii_new * sin(angles)))
    } else {
        radii_new <- round(max(cent) - cent, 8)
        angles <- apply(x, 1, function(y) atan2(y[2], y[1]))
        return(cbind(radii_new * cos(angles), radii_new * sin(angles)))
    }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' constrained stress layout
#'
#' @name layout_constrained_stress
#' @description force-directed graph layout based on stress majorization with variable constrained
#' @param g igraph object
#' @param coord numeric vector. fixed coordinates for dimension specified in `fixdim`.
#' @param fixdim string. which dimension should be fixed. Either "x" or "y".
#' @param weights possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @param mds should an MDS layout be used as initial layout (default: TRUE)
#' @param bbox constrain dimension of output. Only relevant to determine the placement of disconnected graphs
#' @details Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight).
#'
#' The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' 'ggraph' natively supports the layout.
#' @seealso [layout_constrained_stress3D]
#' @return matrix of xy coordinates
#' @references Gansner, E. R., Koren, Y., & North, S. (2004). Graph drawing by stress majorization. *In International Symposium on Graph Drawing* (pp. 239-250). Springer, Berlin, Heidelberg.
#' @export
layout_with_constrained_stress <- function(g, coord, fixdim = "x", weights = NA,
                                           iter = 500, tol = 0.0001, mds = TRUE, bbox = 30) {
    ensure_igraph(g)

    oldseed <- get_seed()
    set.seed(42) # stress is deterministic and produces same result up to translation. This keeps the layout fixed
    on.exit(restore_seed(oldseed))

    fixdim <- match.arg(fixdim, c("x", "y"))
    fixdim <- ifelse(fixdim == "x", 1, 2)

    if (missing(coord)) {
        stop('"coord" is missing with no default.')
    }
    comps <- igraph::components(g, "weak")
    if (comps$no == 1) {
        if (igraph::vcount(g) == 1) {
            return(matrix(c(0, 0), 1, 2))
        } else {
            D <- igraph::distances(g, weights = weights)
            W <- 1 / D^2
            diag(W) <- 0
            n <- igraph::vcount(g)
            xinit <- .init_layout(g, D, mds, n, dim = 2)
            xinit[, fixdim] <- coord
            return(constrained_stress_major(xinit, fixdim, W, D, iter, tol))
        }
    } else {
        return(.component_layouter(
            g = g, weights = weights, comps = comps, dim = 2, mds = mds,
            bbox = bbox, iter = iter, tol = tol, FUN = constrained_stress_major, fixdim = fixdim, coord = coord
        ))
    }
}

#' constrained stress layout in 3D
#'
#' @name layout_constrained_stress3D
#' @description force-directed graph layout based on stress majorization with variable constrained in 3D
#' @param g igraph object
#' @param coord numeric vector. fixed coordinates for dimension specified in `fixdim`.
#' @param fixdim string. which dimension should be fixed. Either "x", "y" or "z".
#' @param weights possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @param mds should an MDS layout be used as initial layout (default: TRUE)
#' @param bbox constrain dimension of output. Only relevant to determine the placement of disconnected graphs
#' @details Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight).
#'
#' This function does not come with direct support for igraph or ggraph.
#'
#' @seealso [layout_constrained_stress]
#' @return matrix of xyz coordinates
#' @references Gansner, E. R., Koren, Y., & North, S. (2004). Graph drawing by stress majorization. *In International Symposium on Graph Drawing* (pp. 239-250). Springer, Berlin, Heidelberg.
#' @export
layout_with_constrained_stress3D <- function(g, coord, fixdim = "x", weights = NA,
                                             iter = 500, tol = 0.0001, mds = TRUE, bbox = 30) {
    ensure_igraph(g)

    oldseed <- get_seed()
    set.seed(42)
    on.exit(restore_seed(oldseed))

    fixdim <- match.arg(fixdim, c("x", "y", "z"))
    fixdim <- ifelse(fixdim == "x", 1, ifelse(fixdim == "y", 2, 3))

    comps <- igraph::components(g, "weak")
    if (comps$no == 1) {
        if (igraph::vcount(g) == 1) {
            return(matrix(c(0, 0, 0), 1, 3))
        } else {
            D <- igraph::distances(g, weights = weights)
            W <- 1 / D^2
            diag(W) <- 0
            n <- igraph::vcount(g)
            xinit <- .init_layout(g, D, mds, n, dim = 3)
            xinit[, fixdim] <- coord
            return(constrained_stress_major3D(xinit, fixdim, W, D, iter, tol))
        }
    } else {
        return(.component_layouter(
            g = g, weights = weights, comps = comps, dim = 3, mds = mds,
            bbox = bbox, iter = iter, tol = tol, FUN = constrained_stress_major3D, fixdim = fixdim, coord = coord
        ))
    }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' Layout with fixed coordinates
#'
#' @name layout_fixed_coords
#' @description force-directed graph layout based on stress majorization with
#' fixed coordinates for some nodes
#' @param g igraph object
#' @param coords numeric n x 2 matrix, where n is the number of nodes. values
#' are either NA or fixed coordinates. coordinates are only calculated for the
#' NA values.
#' @param weights possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @param mds should an MDS layout be used as initial layout (default: TRUE)
#' @param bbox constrain dimension of output. Only relevant to determine the placement of disconnected graphs
#' @details Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight).
#'
#' The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' 'ggraph' natively supports the layout.
#' @seealso [layout_constrained_stress]
#' @return matrix of xy coordinates
#' @examples
#' library(igraph)
#' set.seed(12)
#' g <- sample_bipartite(10, 5, "gnp", 0.5)
#' fxy <- cbind(c(rep(0, 10), rep(1, 5)), NA)
#' xy <- layout_with_fixed_coords(g, fxy)
#' @export
layout_with_fixed_coords <- function(g, coords, weights = NA,
                                     iter = 500, tol = 0.0001, mds = TRUE, bbox = 30) {
    ensure_igraph(g)

    oldseed <- get_seed()
    set.seed(42) # stress is deterministic and produces same result up to translation. This keeps the layout fixed
    on.exit(restore_seed(oldseed))

    if (missing(coords)) {
        stop('"coords" is missing with no default.')
    }
    if (nrow(coords) != igraph::vcount(g) || ncol(coords) != 2) {
        stop("coords has the wrong dimensions")
    }
    if (all(!is.na(coords))) {
        warning("all coordinates fixed")
        return(coords)
    }
    comps <- igraph::components(g, "weak")
    if (comps$no == 1) {
        if (igraph::vcount(g) == 1) {
            return(matrix(c(0, 0), 1, 2))
        } else {
            D <- igraph::distances(g, weights = weights)
            W <- 1 / D^2
            diag(W) <- 0
            n <- igraph::vcount(g)
            xinit <- .init_layout(g, D, mds, n, dim = 2)
            not_na_idx <- which(!is.na(coords), arr.ind = TRUE)
            xinit[not_na_idx] <- coords[not_na_idx]
            return(fixed_stress_major(xinit, coords, W, D, iter, tol))
        }
    } else {
        # return(.component_layouter(
        #     g = g, weights = weights, comps = comps, dim = 2, mds = mds,
        #     bbox = bbox, iter = iter, tol = tol, FUN = constrained_stress_major, fixdim = fixdim, coord = coord
        # ))
        stop("only connected graphs are supported")
    }
}

#' radial focus group layout
#'
#' @description arrange nodes in concentric circles around a focal node according to their distance from the focus and keep predefined groups in the same angle range.
#'
#' @name layout_focus_group
#' @param g igraph object
#' @param v id of focal node to be placed in the center
#' @param group vector indicating grouping of nodes
#' @param shrink shrink the reserved angle range for a group to increase the gaps between groups
#' @param weights possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @details Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight).
#' @seealso [layout_focus]
#' The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' @return matrix of xy coordinates
#' @examples
#' library(igraph)
#' g <- sample_islands(4, 5, 0.8, 2)
#' grp <- as.character(rep(1:4, each = 5))
#' layout_with_focus_group(g, v = 1, group = grp, shrink = 10)
#' @export
layout_with_focus_group <- function(g, v, group, shrink = 10, weights = NA, iter = 500, tol = 0.0001) {
    ensure_igraph(g)
    ensure_connected(g)

    n_grp <- length(unique(group))
    xy <- layout_with_focus(g, v)$xy
    ints <- seq(0, 360, length.out = n_grp + 1)

    for (i in seq_len(n_grp)) {
        xy[group == i, ] <- map_to_angle_range(xy[group == i, ], c(ints[i] + shrink, ints[i + 1] - shrink))
    }
    return(xy)
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' radial centrality group layout
#'
#' @description arranges nodes in concentric circles according to a centrality index and keeping groups within a angle range
#'
#' @name layout_centrality_group
#' @param g igraph object
#' @param cent centrality scores
#' @param group vector indicating grouping of nodes
#' @param shrink shrink the reserved angle range for a group to increase the gaps between groups
#' @param ... additional arguments to layout_with_centrality
#' The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' 'ggraph' natively supports the layout.
#' @seealso [layout_centrality]
#' @return matrix of xy coordinates
#' @examples
#' library(igraph)
#' @export
layout_with_centrality_group <- function(g, cent, group, shrink = 10, ...) {
    ensure_igraph(g)
    ensure_connected(g)
    if (missing(group)) {
        stop('argument "group" is missing with no default.')
    }
    if (missing(cent)) {
        stop('argument "group" is missing with no default.')
    }
    n_grp <- length(unique(group))
    xy <- layout_with_centrality(g, cent, ...)
    ints <- seq(0, 360, length.out = n_grp + 1)

    for (i in seq_len(n_grp)) {
        xy[group == i, ] <- map_to_angle_range(xy[group == i, ], c(ints[i] + shrink, ints[i + 1] - shrink))
    }
    return(xy)
}

#' @useDynLib graphlayouts, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
