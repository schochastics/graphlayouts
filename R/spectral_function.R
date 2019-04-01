#' spectral graph layouts
#' @description Using a set of eigenvectors of matrices associated with a graph as coordinates
#' @name spectral_layout
#' @param g igraph object
#' @param type matrix to be used for spectral decomposition. either 'adjacency' or 'laplacian'
#' @param ev eigenvectors to be used. Either 'smallest' or 'largest'.
#' @details the layout_igraph_* function should not be used directly. It is only used as an argument for 'ggraph'.
#' @return coordinates to be used layouting a graph
#' @examples
#' library(igraph)
#' library(ggraph)
#'
#' g <- sample_gnp(50,0.2)
#'
#' ggraph(g,layout="eigen",type="adjacency",ev="largest")+
#'   geom_edge_link(n=2,edge_colour="grey66")+
#'   geom_node_point(shape=21,fill="grey25",size=5)+
#'   theme_graph()
#'
#' ggraph(g,layout="eigen",type="adjacency",ev="smallest")+
#'   geom_edge_link(n=2,edge_colour="grey66")+
#'   geom_node_point(shape=21,fill="grey25",size=5)+
#'   theme_graph()
#'
#'
#' ggraph(g,layout="eigen",type="laplacian",ev="largest")+
#'   geom_edge_link(n=2,edge_colour="grey66")+
#'   geom_node_point(shape=21,fill="grey25",size=5)+
#'   theme_graph()
#'
#' ggraph(g,layout="eigen",type="laplacian",ev="smallest")+
#'   geom_edge_link(n=2,edge_colour="grey66")+
#'   geom_node_point(shape=21,fill="grey25",size=5)+
#'   theme_graph()
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
