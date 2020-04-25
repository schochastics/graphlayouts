#' @rdname layout_stress
#' @param circular not used
#' @export
layout_igraph_stress <- function(g,weights=NA,iter=500,tol=0.0001,mds=TRUE,bbox=30,circular){
  xy <- layout_with_stress(g,weights,iter,tol,mds,bbox)

  nodes <- data.frame(x=xy[,1],y=xy[,2])
  nodes$circular <- FALSE
  extraData <- as.data.frame(igraph::vertex_attr(g))
  nodes <- cbind(nodes, extraData[, !names(extraData) %in% names(nodes), drop = FALSE])
  nodes
}

#' @rdname layout_focus
#' @param circular not used
#' @export
layout_igraph_focus <- function(g,v,weights=NA,iter=500,tol=0.0001,circular){
  xy <- layout_with_focus(g,v,weights,iter,tol)$xy
  nodes <- data.frame(x=xy[,1],y=xy[,2])
  nodes$circular <- FALSE
  extraData <- as.data.frame(igraph::vertex_attr(g))
  nodes <- cbind(nodes, extraData[, !names(extraData) %in% names(nodes), drop = FALSE])
  nodes
}

#' @rdname layout_centrality
#' @param circular not used
#' @export
layout_igraph_centrality <- function(g,cent,scale=TRUE,iter=500,tol=0.0001,tseq=seq(0,1,0.2),circular){
  xy <- layout_with_centrality(g,cent,scale,iter,tol,tseq)
  nodes <- data.frame(x=xy[,1],y=xy[,2])
  nodes$circular <- FALSE
  extraData <- as.data.frame(igraph::vertex_attr(g))
  nodes <- cbind(nodes, extraData[, !names(extraData) %in% names(nodes), drop = FALSE])
  nodes
}

#' @rdname layout_backbone
#' @param circular not used
#' @export
layout_igraph_backbone <- function(g,keep=0.2,backbone=TRUE,circular){
  xy <- layout_as_backbone(g,keep,backbone)$xy
  nodes <- data.frame(x=xy[,1],y=xy[,2])
  nodes$circular <- FALSE
  extraData <- as.data.frame(igraph::vertex_attr(g))
  nodes <- cbind(nodes, extraData[, !names(extraData) %in% names(nodes), drop = FALSE])
  nodes
}

#' @rdname layout_spectral
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

#' @rdname layout_pmds
#' @param circular not used
#' @export
layout_igraph_pmds <- function(g,pivots,weights=NA,D=NULL,circular){
  xy <- layout_with_pmds(g,pivots,weights)
  nodes <- data.frame(x=xy[,1],y=xy[,2])
  nodes$circular <- FALSE
  extraData <- as.data.frame(igraph::vertex_attr(g))
  nodes <- cbind(nodes, extraData[, !names(extraData) %in% names(nodes), drop = FALSE])
  nodes
}

#' @rdname layout_sparse_stress
#' @param circular not used
#' @export
layout_igraph_sparse_stress <- function(g,pivots,weights=NA,iter=500,circular){
  xy <- layout_with_sparse_stress(g,pivots,weights,iter)
  nodes <- data.frame(x=xy[,1],y=xy[,2])
  nodes$circular <- FALSE
  extraData <- as.data.frame(igraph::vertex_attr(g))
  nodes <- cbind(nodes, extraData[, !names(extraData) %in% names(nodes), drop = FALSE])
  nodes
}

#' @rdname layout_constrained_stress
#' @param circular not used
#' @export
layout_igraph_constrained_stress <- function(g,coord,fixdim="x",weights = NA,
                                             iter = 500,tol = 0.0001,mds = TRUE,bbox = 30,circular){

  xy <- layout_with_constrained_stress(g,coord,fixdim,weights,iter,tol,mds,bbox)

  nodes <- data.frame(x=xy[,1],y=xy[,2])
  nodes$circular <- FALSE
  extraData <- as.data.frame(igraph::vertex_attr(g))
  nodes <- cbind(nodes, extraData[, !names(extraData) %in% names(nodes), drop = FALSE])
  nodes
}

#' @rdname layout_multilevel
#' @param circular not used
#' @export
layout_igraph_multilevel <- function(g, type = "all", FUN1, FUN2,
                                 params1 = NULL, params2 = NULL,
                                 ignore_iso = TRUE,alpha = 35,beta = 45,circular){

  xy <- layout_as_multilevel(g,type,FUN1,FUN2,params1,params2,ignore_iso,alpha,beta)
  nodes <- data.frame(x=xy[,1],y=xy[,2])
  nodes$circular <- FALSE
  extraData <- as.data.frame(igraph::vertex_attr(g))
  nodes <- cbind(nodes, extraData[, !names(extraData) %in% names(nodes), drop = FALSE])
  nodes


}
