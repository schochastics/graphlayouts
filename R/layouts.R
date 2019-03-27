#' Quick plot network
#' @description qgraph is a shortcut designed to obtain a quick view of a network using a stress based layout.
#' @param g igraph object
#' @export
#'
qgraph <- function(g){
  if(!requireNamespace("ggraph", quietly = TRUE)){
    stop("ggraph needed for this function to work. Please install it.", call. = FALSE)
  }
  ggraph::ggraph(g,layout="stress")+
    ggraph::geom_edge_link(edge_colour="grey66")+
    ggraph::geom_node_point(col="black",fill="grey25",shape=21,size=5)+
    ggraph::theme_graph()
}

#' @rdname stress_layout
#' @param circular not used
#' @export
layout_igraph_stress <- function(g,iter=500,tol=0.0001,mds=TRUE,bbox=50,circular){
  xy <- layout_with_stress(g,iter,tol,mds,bbox)

  nodes <- data.frame(x=xy[,1],y=xy[,2])
  nodes$circular <- FALSE
  extraData <- as.data.frame(igraph::vertex_attr(g))
  nodes <- cbind(nodes, extraData[, !names(extraData) %in% names(nodes), drop = FALSE])
  nodes
}

#' @rdname focal_layout
#' @param circular not used
#' @export
layout_igraph_focus <- function(g,v,iter=500,tol=0.0001,circular){
  xy <- layout_with_focus(g,v,iter,tol)
  nodes <- data.frame(x=xy[,1],y=xy[,2])
  nodes$circular <- FALSE
  extraData <- as.data.frame(igraph::vertex_attr(g))
  nodes <- cbind(nodes, extraData[, !names(extraData) %in% names(nodes), drop = FALSE])
  nodes
}

#' @rdname centrality_layout
#' @param circular not used
#' @export
layout_igraph_centrality <- function(g,cent,scale=T,iter=500,tol=0.0001,tseq=seq(0,1,0.2),circular){
  xy <- layout_with_centrality(g,cent,scale,iter,tol,tseq)
  nodes <- data.frame(x=xy[,1],y=xy[,2])
  nodes$circular <- FALSE
  extraData <- as.data.frame(igraph::vertex_attr(g))
  nodes <- cbind(nodes, extraData[, !names(extraData) %in% names(nodes), drop = FALSE])
  nodes
}

#' @rdname backbone_layout
#' @param circular not used
#' @export
layout_igraph_backbone <- function(g,keep=0.2,backbone=T,circular){
  xy <- layout_as_backbone(g,keep,backbone)$xy
  nodes <- data.frame(x=xy[,1],y=xy[,2])
  nodes$circular <- FALSE
  extraData <- as.data.frame(igraph::vertex_attr(g))
  nodes <- cbind(nodes, extraData[, !names(extraData) %in% names(nodes), drop = FALSE])
  nodes
}

#' @rdname spectral_layout
#' @param circular not used
#' @export
layout_igraph_eigen <- function(g,type="laplacian",ev="smallest",circular){
  xy <- layout_with_eigen(g,type,ev)
  nodes <- data.frame(x=xy[,1],y=xy[,2])
  nodes$circular <- FALSE
  extraData <- as.data.frame(igraph::vertex_attr(g))
  nodes <- cbind(nodes, extraData[, !names(extraData) %in% names(nodes), drop = FALSE])
  nodes
}
