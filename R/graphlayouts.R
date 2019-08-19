#' graphlayouts: layout algorithms for network visualizations
#'
#' @description The package implements several new layout algorithms to visualize networks.
#' Most are based on the concept of stress majorization. Some more specific algorithms allow to emphasize
#' hidden group structures in networks or focus on specific nodes. The package is best used in conjunction with
#' ggraph.
#'
#' @details Some features of the package are:
#'
#' \itemize{
#' \item `layout_with_stress()` is a state of the art deterministic layout algorithms.
#' \item `layout_as_backbone()` uncovers hidden group structures (if they exist) by emphasizing strongly embedded edges.
#' \item `layout_with_focus()` and `layout_with_centrality()` produce concentric layouts with a focal or most central nodes in the center.
#' \item `layout_with_eigen()` implements some layout algorithms on the basis of eigenvectors
#' \item `layout_with_sparse_stress()` sparse stress for large graphs
#' \item `layout_with_pmds()` pivot MDS for large graphs.
#' \item `layout_as_dynamic()` for longitudinal network data
#' }
#'
#' A detailed tutorial can be found \href{http://mr.schochastics.net/netVizR.html}{here}.
#' @docType package
#' @name graphlayouts
NULL
