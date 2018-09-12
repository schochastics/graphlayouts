#' Stress majorization graph layout
#'
#' @param g igraph object
#' @param dim dimension of output coordinates. default is 2
#' @param iter number of iterations
#' @param tol stoping criterion
#'
#' @return coordinates to be used layouting a graph
#' @export
#'
stress_majorization <- function(g,dim=2,iter=500,tol=0.0001){
  D <- igraph::distances(g)
  W <- 1/D^2
  diag(W) <- 0
  xinit <- matrix(stats::runif(igraph::vcount(g)*dim,0,1),igraph::vcount(g),dim)
  x <- stress_major(xinit,W,D,dim,iter,tol)
  x
}

#' @useDynLib smglr
#' @importFrom Rcpp sourceCpp
NULL
