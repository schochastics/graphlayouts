#' dynamic graph layout
#' @description Create layouts for longitudinal networks.
#' @name layout_dynamic
#' @param gList list of igraph objects. Each network must contain the same set of nodes.
#' @param weights possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @param alpha weighting of reference layout. See details.
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @details The reference layout is calculated based on the union of all graphs. The parameter alpha controls the influence of the reference layout.
#' For alpha=1, only the reference layout is used and all graphs have the same layout. For alpha=0, the stress layout of each individual graph is used. Values in-between interpolate between the two layouts.
#'
#' Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight).
#' @return list of coordinates for each graph
#' @references Brandes, U. and Indlekofer, N. and Mader, M. (2012). Visualization methods for longitudinal social networks and stochastic actor-oriented modeling. *Social Networks* 34 (3) 291-308
#' @examples
#' library(igraph)
#' g1 <- sample_gnp(20,0.2)
#' g2 <- sample_gnp(20,0.2)
#' g3 <- sample_gnp(20,0.2)
#'
#' xy <- layout_as_dynamic(list(g1,g2,g3))
#'
#' # layout for first network
#' xy[[1]]
#' @export
layout_as_dynamic <- function(gList,weights = NA, alpha = 0.5,iter = 500,tol = 1e-04){
  check_networks <- vapply(gList,FUN = function(x) igraph::is_igraph(x),FUN.VALUE = FALSE)
  if(!all(check_networks)){
    stop("'gList' must be a list of igraph objects.")
  }
  #prepare reference layout
  g <- Reduce("%u%",gList)
  check_nodes <- vapply(gList,FUN = function(x) igraph::vcount(x)==igraph::vcount(g),FUN.VALUE = FALSE)
  if(!all(check_nodes)){
    stop("all nodes must be present in each network")
  }
  n <- igraph::vcount(g)
  DList <- lapply(gList,igraph::distances,weights = weights)
  DList <- adjust_dist(DList)
  Dmean <- Reduce('+', DList)/length(DList)
  Dvar <- Reduce('+',lapply(DList, function(x) (x-Dmean)^2))/length(DList)
  W <- 1/Dmean^2+1/(1+Dvar)
  diag(W) <- 0

  #calculate reference layout
  rmat <- matrix(stats::runif(n*2,-0.1,0.1),n,2)
  xinit <- igraph::layout_with_mds(g) + rmat
  xref <- stress_major(xinit,W,Dmean,iter,tol)

  xycoords <-vector("list",length(gList))
  for(i in 1:length(gList)){
    D <- DList[[i]]
    W <- 1/D^2
    diag(W) <- 0
    if(i==1){
      xycoords[[i]] <- stress_major(xref,W,D,iter,tol)
    } else{
      xycoords[[i]] <- stress_major(xycoords[[i-1]],W,D,iter,tol)
    }
    xycoords[[i]] <- (1-alpha)*xycoords[[i]]+alpha*xref
  }
  xycoords
}

adjust_dist <- function(DList){
  n <- nrow(DList[[1]])
  for(i in 1:n){
    for(j in 1:n){
      for(k in 1:length(DList)){
        if(is.infinite(DList[[k]][i,j])){
          lastD <- Inf
          for(l in seq((k-1),1)){
            if(l==0){
              next()
            }
            if(!is.infinite(DList[[l]][i,j])){
              lastD <- DList[[l]][i,j]
              tlast <- l
              break()
            }
          }
          nextD <- Inf
          for(l in seq((k+1),length(DList))){
            if(l>length(DList)){
              break()
            }
            if(!is.infinite(DList[[l]][i,j])){
              nextD <- DList[[l]][i,j]
              tnext <- l
              break()
            }
          }
          if(!is.infinite(lastD) & !is.infinite(nextD)){
            beta <- (k-tlast)/(tnext-tlast)
            DList[[k]][i,j] <- (1-beta)*lastD+beta*nextD+1
          } else if(is.infinite(lastD) & !is.infinite(nextD)){
            DList[[k]][i,j] <- nextD+1
          } else if(!is.infinite(lastD) & is.infinite(nextD)){
            DList[[k]][i,j] <- lastD+1
          } else{
            DList[[k]][i,j] <- sqrt(n)
          }
        }
      }
    }
  }
  DList
}
