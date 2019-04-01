#' @title Manipulate graph
#' @description functions to manipulate graphs
#'
#' @param g igraph object
#' @param attr edge attribute name used to sort edges
#' @param desc logical. sort in descending (default) or ascending order
#' @details `reorder_edges()` allows to reorder edges according to an attribute so that edges are
#' drawn in the given order.
#' @name graph_manipulate
#' @return manipulated graphs
#' @examples
#' library(igraph)
#' library(ggraph)
#'
#' g <- sample_gnp(10,0.5)
#' E(g)$attr <- 1:ecount(g)

#' ggraph(g,layout="stress")+
#'   geom_edge_link(aes(col=attr))+
#'   geom_node_point()+
#'   theme_graph()+
#'   theme(legend.position="none")
#'
#' gn <- reorder_edges(g,"attr")
#' ggraph(gn,layout="stress")+
#'   geom_edge_link(aes(col=attr))+
#'   geom_node_point()+
#'   theme_graph()+
#'   theme(legend.position="none")
#'
#' @author David Schoch
NULL

#' @rdname graph_manipulate
#' @export

reorder_edges <- function(g,attr,desc=TRUE){
  if(!"name"%in%igraph::vertex_attr_names(g)){
    igraph::V(g)$name <- 1:igraph::vcount(g)
  }
  edges_df <- igraph::as_data_frame(g,what="edges")
  edges_df <- edges_df[order(edges_df[[attr]],decreasing = desc),]

  vertices <- igraph::as_data_frame(g,what="vertices")
  vattrs <- igraph::vertex_attr_names(g)
  idname <- which(vattrs=="name")
  vertices <- vertices[,c(idname,setdiff(1:length(vattrs),idname))]

  gn <- igraph::graph_from_data_frame(d = edges_df,
                                      directed = igraph::is.directed(g),
                                      vertices = vertices)
  gn
}

