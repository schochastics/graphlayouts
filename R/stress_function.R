#' Stress majorization graph layout
#'
#' @param g igraph object
#' @param iter number of iterations
#' @param tol stoping criterion
#' @param mds should an MDS layout be used as initial layout (default: TRUE)
#'
#' @return coordinates to be used layouting a graph
#' @export
#'
stress_majorization <- function(g,iter=500,tol=0.0001,mds=TRUE){
  D <- igraph::distances(g,weights = NA)
  W <- 1/D^2
  diag(W) <- 0
  n <- igraph::vcount(g)
  if(!mds){
    xinit <- matrix(stats::runif(n*2,0,1),igraph::vcount(g),2)
  } else{
    rmat <- matrix(stats::runif(n*2,-0.1,0.1),igraph::vcount(g),2)
    xinit <- igraph::layout_with_mds(g) + rmat
  }
  x <- stress_major(xinit,W,D,iter,tol)
  x
}

#' @useDynLib smglr
#' @importFrom Rcpp sourceCpp
NULL
