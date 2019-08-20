#' spectral graph layouts
#' @description Using a set of eigenvectors of matrices associated with a graph as coordinates
#' @name layout_spectral
#' @param g igraph object
#' @param type matrix to be used for spectral decomposition. either 'adjacency' or 'laplacian'
#' @param ev eigenvectors to be used. Either 'smallest' or 'largest'.
#' @details The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' 'ggraph' natively supports the layout.
#' @author David Schoch
#' @return matrix of xy coordinates
#' @examples
#' library(igraph)
#'
#' g <- sample_gnp(50,0.2)
#'
#' xy <- layout_with_eigen(g,type = "adjacency",ev = "largest")
#'
#' xy <- layout_with_eigen(g,type = "adjacency",ev = "smallest")
#'
#' xy <- layout_with_eigen(g,type = "laplacian",ev = "largest")
#'
#' xy <- layout_with_eigen(g,type = "laplacian",ev = "smallest")
#' @export

layout_with_eigen <- function(g,type="laplacian",ev="smallest"){
  if(!igraph::is_connected(g)){
    stop("g must be connected")
  }
  if(igraph::is_directed(g)){
    warning("g is directed. undirected version is used for the layout.")
    g <- igraph::as.undirected(g)
  }
  if(!type%in%c("laplacian","adjacency")){
    stop("type must be one of 'laplacian' or 'adjacency'")
  }
  if(!ev%in%c("largest","smallest")){
    stop("ev must be one of 'smallest' or 'largest'")
  }
  n <- igraph::vcount(g)
  if(type=="adjacency"){
    A <- igraph::get.adjacency(g,type="both")
  } else{
    A <- igraph::laplacian_matrix(g)
  }
  sA <- eigen(A)
  if(ev=="largest"){
    xy <- sA$vectors[,1:2]
  } else if(ev=="smallest" & type=="adjacency"){
    xy <- sA$vectors[,(n-1):n]
  } else{
    xy <- sA$vectors[,(n-2):(n-1)]
  }
  xy
}
