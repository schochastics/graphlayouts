#' Stress majorization graph layout
#'
#' @param g igraph object
#' @param iter number of iterations
#' @param tol stoping criterion
#' @param mds should an MDS layout be used as initial layout (default: TRUE)
#' @param bbox constrain dimension of output
#'
#' @return coordinates to be used layouting a graph
#' @export
#'
stress_majorization <- function(g,iter=500,tol=0.0001,mds=TRUE,bbox=50){

  comps <- igraph::components(g,"weak")
  if(comps$no>1){
    lg <- list()

    for (i in 1:comps$no){
      sg <- igraph::induced_subgraph(g,comps$membership==i)
      n <- igraph::vcount(sg)
      if(n==1){
        lg[[i]] <- matrix(c(0,0),1,2,byrow = T)
        next()
      }
      if(n==2){
        lg[[i]] <- matrix(c(0,0,1,0),2,2,byrow = T)
        next()
      }

      D <- igraph::distances(sg,weights = NA)
      W <- 1/D^2
      diag(W) <- 0
      if(!mds){
        xinit <- matrix(stats::runif(n*2,0,1),n,2)
      } else{
        rmat <- matrix(stats::runif(n*2,-0.1,0.1),n,2)
        xinit <- igraph::layout_with_mds(sg) + rmat
      }
      lg[[i]] <- stress_major(xinit,W,D,iter,tol)
    }
    lg <- lapply(lg,mv_to_null)
    # lg <- lapply(lg,function(x) x/max(x[,1]))

    ldhw <- do.call("rbind",lapply(lg,get_bbox))[,3:4]+1

    p <- order(comps$csize,decreasing=T)-1
    m <- c(bbox,bbox)

    l2d <- gbp::gbp2d_solver_dpp(p, t(ldhw), m)
    for(i in 1:comps$no){
      lg[[i]][,1] <- lg[[i]][,1]+l2d$it[1,i]
      lg[[i]][,2] <- lg[[i]][,2]+l2d$it[2,i]
    }
    x <- do.call("rbind",lg)

  } else{
    if(igraph::vcount(g)==1){
      x <- matrix(c(0,0),1,2)
    } else{
      D <- igraph::distances(g,weights = NA)
      W <- 1/D^2
      diag(W) <- 0
      n <- igraph::vcount(g)
      if(!mds){
        xinit <- matrix(stats::runif(n*2,0,1),n,2)
      } else{
        rmat <- matrix(stats::runif(n*2,-0.1,0.1),n,2)
        xinit <- igraph::layout_with_mds(g) + rmat
      }
      x <- stress_major(xinit,W,D,iter,tol)
    }
  }
  x
}


#-------------------------------------------------------------------------------
# helper functions
#-------------------------------------------------------------------------------

get_bbox <- function(xy){
  lbottom <- c(min(xy[,1]),min(xy[,2]))
  rtop <- c(max(xy[,1]),max(xy[,2]))
  c(lbottom,rtop)
}

mv_to_null <- function(xy){
  bbox <- get_bbox(xy)
  xy[,1] <- xy[,1]-bbox[1]
  xy[,2] <- xy[,2]-bbox[2]
  xy
}


#' @useDynLib smglr
#' @importFrom Rcpp sourceCpp
NULL
