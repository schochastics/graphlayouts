#' pivot MDS graph layout
#'
#' @name layout_pmds
#' @description similar to \link[igraph]{layout_with_mds} but uses only a small set of pivots for MDS. Considerably faster than MDS and thus applicable for larger graphs.
#' @param g igraph object
#' @param pivots number of pivots
#' @param D precomputed distances from pivots to all nodes (if available, default: NULL)
#' @param dim dimensionality of layout (defaults to 2)
#' @param weights possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @details Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight)
#'
#' The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' 'ggraph' natively supports the layout.
#' @author David Schoch
#' @return matrix of coordinates
#' @references Brandes, U. and Pich, C. (2006). Eigensolver Methods for Progressive Multidimensional Scaling of Large Data. In *International Symposium on Graph Drawing* (pp. 42-53). Springer
#' @examples
#' \dontrun{
#' library(igraph)
#' library(ggraph)
#'
#' g <- sample_gnp(1000, 0.01)
#'
#' xy <- layout_with_pmds(g, pivots = 100)
#' }
#' @export
layout_with_pmds <- function(g, pivots, weights = NA, D = NULL, dim = 2) {
    ensure_igraph(g)
    ensure_connected(g)

    if (missing(pivots) && is.null(D)) {
        stop('argument "pivots" is missing, with no default.')
    }
    if (!missing(pivots)) {
        if (pivots > igraph::vcount(g)) {
            stop('"pivots" must be less than the number of nodes in the graph.')
        }
    }
    if (is.null(D)) {
        pivs <- sample(1:igraph::vcount(g), pivots)
        D <- t(igraph::distances(g, v = pivs, weights = weights))
    }
    cmean <- colMeans(D^2)
    rmean <- rowMeans(D^2)
    Dmat <- D^2 - outer(rmean, cmean, function(x, y) x + y) + mean(D^2)
    sl2 <- svd(Dmat)

    xy <- (Dmat %*% sl2$v[, 1:dim])
    xy
}


#' sparse stress graph layout
#'
#' @name layout_sparse_stress
#' @description stress majorization for larger graphs based on a set of pivot nodes.
#' @param g igraph object
#' @param pivots number of pivots
#' @param weights ignored
#' @param iter number of iterations during stress optimization
#' @details The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' 'ggraph' natively supports the layout.
#' @author David Schoch
#' @return matrix of xy coordinates
#' @references Ortmann, M. and Klimenta, M. and Brandes, U. (2016). A Sparse Stress Model. https://arxiv.org/pdf/1608.08909.pdf
#' @examples
#' \dontrun{
#' library(igraph)
#' library(ggraph)
#'
#' g <- sample_gnp(1000, 0.005)
#'
#' ggraph(g, layout = "sparse_stress", pivots = 100) +
#'     geom_edge_link0(edge_colour = "grey66") +
#'     geom_node_point(shape = 21, fill = "grey25", size = 5) +
#'     theme_graph()
#' }
#' @export

layout_with_sparse_stress <- function(g, pivots, weights = NA, iter = 500) {
    ensure_igraph(g)
    ensure_connected(g)
    if (!all(is.na(weights))) {
        warning("weights are not supported. unweighted graph is used instead.")
    }
    if (missing(pivots)) {
        stop('argument "pivots" is missing, with no default.')
    }
    if (pivots > igraph::vcount(g)) {
        stop('"pivots" must be less than the number of nodes in the graph.')
    }
    comps <- igraph::components(g)
    if (comps$no == 1) {
        prep <- .sparse_prepare(g, pivots)
        A <- igraph::as_adj(g, type = "both", sparse = TRUE)
        return(sparseStress(prep$y, prep$D, prep$RpL, prep$pivs, A, iter))
    } else {
        # TBD
    }
}

.sparse_prepare <- function(g, pivots) {
    pivs <- sample(1:igraph::vcount(g), pivots)

    D <- t(igraph::distances(g, v = pivs, weights = NA))
    Rp <- apply(D, 1, which.min)
    y <- layout_with_pmds(g, pivots, D = D, weights = NA)

    # rescale
    el <- igraph::as_edgelist(g, names = FALSE)
    norm1 <- sum(sqrt((y[el[, 1], 1] - y[el[, 2], 1])^2 + (y[el[, 1], 2] - y[el[, 2], 2])^2))

    y <- y * (igraph::ecount(g) / norm1)

    RpL <- lapply(seq_along(pivs), function(x) which(Rp == x) - 1)
    pivs <- pivs - 1
    list(RpL = RpL, pivs = pivs, y = y, D = D)
}
