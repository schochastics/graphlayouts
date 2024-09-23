#' @title Metro Map Layout
#' @description Metro map layout based on multicriteria optimization
#' @param object original graph
#' @param xy initial layout of the original graph
#' @param l desired multiple of grid point spacing. (l*gr determines desired edge length)
#' @param gr grid spacing. (l*gr determines desired edge length)
#' @param w weight vector for criteria (see details)
#' @param bsize number of grid points a station can move away rom its original position
#' @details The function optimizes the following five criteria using a hill climbing algorithm:
#' - *Angular Resolution Criterion*: The angles of incident edges at each station should be maximized, because if there is only a small angle between any two adjacent edges, then it can become difficult to distinguish between them
#' - *Edge Length Criterion*: The edge lengths across the whole map should be approximately equal to ensure regular spacing between stations. It is based on the preferred multiple, l, of the grid spacing, g. The purpose of the criterion is to penalize edges that are longer than or shorter than lg.
#' - *Balanced Edge Length Criterion*: The length of edges incident to a particular station should be similar
#' - *Line Straightness Criterion*: (not yet implemented) Edges that form part of a line should, where possible, be co-linear either side of each station that the line passes through
#' - *Octilinearity Criterion*: Each edge should be drawn horizontally, vertically, or diagonally at 45 degree, so we penalize edges that are not at a desired angle
#' @return new coordinates for stations
#' @references
#' Stott, Jonathan, et al. "Automatic metro map layout using multicriteria optimization." IEEE Transactions on Visualization and Computer Graphics 17.1 (2010): 101-114.
#' @author David Schoch
#' @examples
#' # the algorithm has problems with parallel edges
#' library(igraph)
#' g <- simplify(metro_berlin)
#' xy <- cbind(V(g)$lon, V(g)$lat) * 100
#'
#' # the algorithm is not very stable. try playing with the parameters
#' \dontrun{
#' xy_new <- layout_as_metromap(g, xy, l = 2, gr = 0.5, w = c(100, 100, 1, 1, 100), bsize = 35)
#' }
#' @export
layout_as_metromap <- function(object, xy, l = 2, gr = 0.0025, w = rep(1, 5), bsize = 5) {
    adj <- as_adj_list1(object)
    adj <- lapply(adj, function(x) x - 1)
    adj_deg2 <- adj[unlist(lapply(adj, length)) == 2]
    el <- igraph::get.edgelist(object, FALSE) - 1

    xy <- snap_to_grid(xy, gr)

    bbox <- station_bbox(xy, bsize, gr)

    xy_new <- layout_as_metro_iter(adj, el, adj_deg2, xy, bbox, l, gr, w, bsize)
    xy_new
}

# helper ----
snap_to_grid <- function(xy, gr) {
    xmin <- min(xy[, 1])
    xmax <- max(xy[, 1])
    ymin <- min(xy[, 2])
    ymax <- max(xy[, 2])

    deltax <- seq(xmin - 4 * gr, xmax + 4 * gr, by = gr)
    deltay <- seq(ymin - 4 * gr, ymax + 4 * gr, by = gr)

    xdiff <- outer(xy[, 1], deltax, function(x, y) abs(x - y))
    ydiff <- outer(xy[, 2], deltay, function(x, y) abs(x - y))

    xy_new <- cbind(deltax[apply(xdiff, 1, which.min)], deltay[apply(ydiff, 1, which.min)])
    dups <- duplicated(xy_new)
    while (any(dups)) {
        xy_new[which(dups), ] <- xy_new[which(dups), ] + c(sample(c(1, -1), 1) * gr, sample(c(1, -1), 1) * gr)
        dups <- duplicated(xy_new)
    }
    xy_new
}

station_bbox <- function(xy, bsize, gr) {
    cbind(xy - bsize * gr, xy + bsize * gr)
}

as_adj_list1 <- function(g) {
    n <- igraph::vcount(g)
    lapply(1:n, function(i) {
        x <- g[[i]][[1]]
        attr(x, "env") <- NULL
        attr(x, "graph") <- NULL
        class(x) <- NULL
        x
    })
}
