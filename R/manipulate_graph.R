#' @title Manipulate graph
#' @description functions to manipulate graphs, such as changing order of edges
#'
#' @param g igraph object
#' @param attr attribute to sort edges
#' @param desc logical. sort in descending (default) or ascending order
#' @name graph_manipulate
#' @return manipulated graphs
#' @author David Schoch
NULL

#' @rdname graph_manipulate
#' @export

reorder_edges <- function(g,attr,desc=T){
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

