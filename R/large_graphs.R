#' pivot MDS graph layout
#'
#' @name pivotMDS
#' @description similar to MDS layout but uses only a small set of pivots for MDS. Considerably faster than MDS and thus applicable for larger graphs.
#' @param g igraph object
#' @param pivots number of pivots
#' @param weights Possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @details the layout_igraph_* function should not be used directly. It is only used as an argument for 'ggraph'.
#' Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight)
#' @return coordinates to be used layouting a graph
#' @references Brandes, U. and Pich, C. (2006). Eigensolver Methods for Progressive Multidimensional Scaling of Large Data. In *International Symposium on Graph Drawing* (pp. 42-53). Springer
#' @examples
#' \dontrun{
#' library(igraph)
#' library(ggraph)
#' g <- sample_gnp(1000,0.01)
#' ggraph(g,layout="pMDS",pivots=250)+
#'   geom_edge_link(n=2,edge_colour="grey66")+
#'   geom_node_point(shape=21,fill="grey25",size=5)+
#'   theme_graph()
#'}
#' @export

layout_with_pmds <- function(g,pivots,weights=NA){
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  if(is.null(pivots)){
    stop('argument "pivots" is missing, with no default')
  }
  pivs <- sample(1:igraph::vcount(g),pivots)
  D <- igraph::distances(g,to=pivs,weights = weights)
  cmean <- colMeans(D^2)
  rmean <- rowMeans(D^2)
  Dmat <- D^2-outer(rmean,cmean, function(x,y) x+y)+mean(D^2)
  sl2 <- svd(Dmat)

  xy <- (Dmat%*%sl2$v[,1:2])
  xy
}


#' sparse stress graph layout
#'
#' @name sparseStress
#' @description stress majorization for larger graphs
#' @param g igraph object
#' @param pivots number of pivots
#' @details the layout_igraph_* function should not be used directly. It is only used as an argument for 'ggraph'. Edge weights are not suported yet.
#' @return coordinates to be used layouting a graph
#' @references Ortmann, M. and Klimenta, M. and Brandes, U. (2016).A Sparse Stress Model. https://arxiv.org/pdf/1608.08909.pdf
#' @examples
#' \dontrun{
#' library(igraph)
#' library(ggraph)
#' g <- sample_gnp(1000,0.005)
#' ggraph(g,layout="sparseStress",pivots=250)+
#'   geom_edge_link(n=2,edge_colour="grey66")+
#'   geom_node_point(shape=21,fill="grey25",size=5)+
#'   theme_graph()
#'}
#' @export

layout_with_sparseStress <- function(g,pivots){
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  pivs <- sample(1:igraph::vcount(g),pivots)

  D <- igraph::distances(g,to=pivs,weights = NA)
  Rp <- apply(D,1,which.min)
  y <- layout_with_pmds(g,pivots)

  #rescale
  el <- igraph::get.edgelist(g,F)
  norm1 <- sum(sqrt((y[el[,1],1]-y[el[,2],1])^2+(y[el[,1],2]-y[el[,2],2])^2))
  n <- igraph::vcount(g)
  y <- y*(igraph::ecount(g)/norm1)+(matrix(runif(n*2,-0.2,0.2),n,2))

  adjL <- igraph::get.adjlist(g,"all")
  RpL <- lapply(1:length(pivs),function(x) which(Rp==x)-1)
  pivs <- pivs-1
  adjL <- lapply(adjL,function(x) x-1)

  xy <- sparseStress(y,D,RpL,pivs,adjL)
  xy
}
