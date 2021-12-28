#' UMAP graph layouts
#' @description Using the UMAP dimensionality reduction algorithm as a graph layout
#' @name layout_umap
#' @param g igraph object
#' @param pivots if not NULL, number of pivot nodes to use for distance calculation (for large graphs).
#' @param ... additional parameters for umap. See the `?uwot::umap` for help.
#' @details The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'. UMAP can be tuned by many different parameters. Refer to the documentation at https://github.com/jlmelville/uwot for help
#' @references
#' McInnes, Leland, John Healy, and James Melville. "Umap: Uniform manifold approximation and projection for dimension reduction." arXiv preprint arXiv:1802.03426 (2018).
#' @author David Schoch
#' @return matrix of xy coordinates
#' @examples
#' library(igraph)
#'
#' g <- sample_islands(10,20,0.6,10)
#' xy <- layout_with_umap(g,min_dist = 0.5)
#' @export

layout_with_umap <- function(g,pivots=NULL,...){
  if(!requireNamespace("uwot", quietly = TRUE)){
    stop("uwot is needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  if(is.null(pivots)){
    D <- igraph::distances(g)
  } else{
    pivs <- sample(1:igraph::vcount(g),pivots)
    D <- t(igraph::distances(g,v=pivs))
  }
  uwot::umap(D,...)
}
