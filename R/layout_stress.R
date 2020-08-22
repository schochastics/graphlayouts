#' stress majorization layout
#'
#' @name layout_stress
#' @rdname layout_stress
#' @description force-directed graph layout based on stress majorization.
#' @param g igraph object
#' @param weights possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @param mds should an MDS layout be used as initial layout (default: TRUE)
#' @param bbox constrain dimension of output. Only relevant to determine the placement of disconnected graphs
#' @details Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight).
#'
#' The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' 'ggraph' natively supports the layout.
#' @seealso [layout_stress3D]
#' @return matrix of xy coordinates
#' @references Gansner, E. R., Koren, Y., & North, S. (2004). Graph drawing by stress majorization. *In International Symposium on Graph Drawing* (pp. 239-250). Springer, Berlin, Heidelberg.
#' @examples
#' library(igraph)
#' library(ggraph)
#' set.seed(665)
#'
#' g <- sample_pa(100,1,1,directed = FALSE)
#'
#' # calculate layout manually
#' xy <- layout_with_stress(g)
#'
#' # use it with ggraph
#' \dontrun{
#' ggraph(g,layout = "stress")+
#'   geom_edge_link0(edge_width = 0.2,colour = "grey")+
#'   geom_node_point(col = "black",size = 0.3)+
#'   theme_graph()
#'  }
#' @export
layout_with_stress <- function(g,weights = NA, iter = 500,tol = 0.0001,mds = TRUE,bbox = 30){
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  if (exists(".Random.seed", .GlobalEnv)){
    oldseed <- .GlobalEnv$.Random.seed
  }
  else{
    oldseed <- NULL
  }
  set.seed(42)  #stress is deterministic and produces same result up to translation. This keeps the layout fixed
  on.exit(restore_seed(oldseed))

  comps <- igraph::components(g,"weak")
  if(comps$no>1){
    lg <- list()
    node_order <- c()
    for (i in 1:comps$no){
      sg <- igraph::induced_subgraph(g,comps$membership==i)
      n <- igraph::vcount(sg)
      node_order <- c(node_order,which(comps$membership==i))
      if(n==1){
        lg[[i]] <- matrix(c(0,0),1,2,byrow = TRUE)
        next()
      }
      if(n==2){
        lg[[i]] <- matrix(c(0,0,1,0),2,2,byrow = TRUE)
        next()
      }

      D <- igraph::distances(sg,weights = weights)
      W <- 1/D^2
      diag(W) <- 0
      if(!mds){
        xinit <- matrix(stats::runif(n*2,0,1),n,2)
      } else{
        rmat <- matrix(stats::runif(n*2,-0.1,0.1),n,2)
        if(igraph::vcount(sg)<=100){
          xinit <- igraph::layout_with_mds(sg) + rmat
        } else{
          xinit <- layout_with_pmds(sg,D = D[,sample(1:igraph::vcount(sg),100)]) + rmat
        }
      }
      lg[[i]] <- stress_major(xinit,W,D,iter,tol)
    }
    lg <- lapply(lg,mv_to_null)
    p <- order(comps$csize)
    curx <- 0
    cury <- 0
    maxy <- 0
    for(comp in p){
      if(curx+max(lg[[comp]][,1])>bbox){
        curx <- 0
        cury <- maxy+1
      }
      lg[[comp]][,1] <- lg[[comp]][,1] + curx
      lg[[comp]][,2] <- lg[[comp]][,2] + cury
      curx <- max(lg[[comp]][,1])+1
      maxy <- max(c(maxy,max(lg[[comp]][,2])))
    }
    x <- do.call("rbind",lg)
    x <- x[order(node_order),]

  } else{
    if(igraph::vcount(g)==1){
      x <- matrix(c(0,0),1,2)
    } else{
      D <- igraph::distances(g,weights = weights)
      W <- 1/D^2
      diag(W) <- 0
      n <- igraph::vcount(g)
      if(!mds){
        xinit <- matrix(stats::runif(n*2,0,1),n,2)
      } else{
        rmat <- matrix(stats::runif(n*2,-0.1,0.1),n,2)
        if(igraph::vcount(g)<=100){
          xinit <- igraph::layout_with_mds(g) + rmat
        } else{
          xinit <- layout_with_pmds(g,D = D[,sample(1:(igraph::vcount(g)),100)]) + rmat
        }

      }
      x <- stress_major(xinit,W,D,iter,tol)
    }
  }
  x
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' stress majorization layout in 3D
#'
#' @name layout_stress3D
#' @rdname layout_stress3D
#' @description force-directed graph layout based on stress majorization in 3D.
#' @param g igraph object
#' @param weights possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @param mds should an MDS layout be used as initial layout (default: TRUE)
#' @param bbox constrain dimension of output. Only relevant to determine the placement of disconnected graphs
#' @details Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight).
#'
#' @return matrix of xyz coordinates
#' @seealso [layout_stress]
#' @references Gansner, E. R., Koren, Y., & North, S. (2004). Graph drawing by stress majorization. *In International Symposium on Graph Drawing* (pp. 239-250). Springer, Berlin, Heidelberg.
#' @export
layout_with_stress3D <- function(g,weights = NA, iter = 500,tol = 0.0001,mds = TRUE,bbox = 30){
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  if (exists(".Random.seed", .GlobalEnv)){
    oldseed <- .GlobalEnv$.Random.seed
  }
  else{
    oldseed <- NULL
  }
  set.seed(42)  #stress is deterministic and produces same result up to translation. This keeps the layout fixed
  on.exit(restore_seed(oldseed))

  comps <- igraph::components(g,"weak")
  if(comps$no>1){
    lg <- list()
    node_order <- c()
    for (i in 1:comps$no){
      sg <- igraph::induced_subgraph(g,comps$membership==i)
      n <- igraph::vcount(sg)
      node_order <- c(node_order,which(comps$membership==i))
      if(n==1){
        lg[[i]] <- matrix(c(0,0,0),1,3,byrow = TRUE)
        next()
      }
      if(n==2){
        lg[[i]] <- matrix(c(0,0,0,1,0,0),2,3,byrow = TRUE)
        next()
      }

      D <- igraph::distances(sg,weights = weights)
      W <- 1/D^2
      diag(W) <- 0
      if(!mds){
        xinit <- matrix(stats::runif(n*3,0,1),n,3)
      } else{
        rmat <- matrix(stats::runif(n*3,-0.1,0.1),n,3)
        if(igraph::vcount(sg)<=100){
          xinit <- igraph::layout_with_mds(sg,dim = 3) + rmat
        } else{
          n <- igraph::vcount(g)
          pivs <- sample(1:n,100)
          D1 <- D[,pivs]
          cmean <- colMeans(D1^2)
          rmean <- rowMeans(D1^2)
          Dmat <- D1^2-outer(rmean,cmean, function(x,y) x+y)+mean(D1^2)
          sl2 <- svd(Dmat)
          rmat <- matrix(stats::runif(n*3,-0.1,0.1),n,3)
          xinit <- (Dmat%*%sl2$v[,1:3]) + rmat
        }
      }
      lg[[i]] <- stress_major3D(xinit,W,D,iter,tol)
    }
    lg <- lapply(lg,mv_to_null)
    p <- order(comps$csize)
    curx <- 0
    cury <- 0
    maxy <- 0
    for(comp in p){
      if(curx+max(lg[[comp]][,1])>bbox){
        curx <- 0
        cury <- maxy+1
      }
      lg[[comp]][,1] <- lg[[comp]][,1] + curx
      lg[[comp]][,2] <- lg[[comp]][,2] + cury
      curx <- max(lg[[comp]][,1])+1
      maxy <- max(c(maxy,max(lg[[comp]][,2])))
    }
    x <- do.call("rbind",lg)
    x <- x[order(node_order),]

  } else{
    if(igraph::vcount(g)==1){
      x <- matrix(c(0,0),1,2)
    } else{
      D <- igraph::distances(g,weights = weights)
      W <- 1/D^2
      diag(W) <- 0
      n <- igraph::vcount(g)
      if(!mds){
        xinit <- matrix(stats::runif(n*2,0,1),n,2)
      } else{
        rmat <- matrix(stats::runif(n*2,-0.1,0.1),n,2)
        if(igraph::vcount(g)<=100){
          xinit <- igraph::layout_with_mds(g) + rmat
        } else{
          xinit <- layout_with_pmds(g,D = D[,sample(1:(igraph::vcount(g)),100)]) + rmat
        }

      }
      x <- stress_major(xinit,W,D,iter,tol)
    }
  }
  x
}


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' radial focus layout
#'
#' @description arrange nodes in concentric circles around a focal node according to their distance from the focus.
#'
#' @name layout_focus
#' @param g igraph object
#' @param v id of focal node to be placed in the center
#' @param weights possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @details Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight).
#'
#' The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' 'ggraph' natively supports the layout.
#' @return a list containing xy coordinates and the distances to the focal node
#' @references Brandes, U., & Pich, C. (2011). More flexible radial layout. *Journal of Graph Algorithms and Applications*, 15(1), 157-173.
#' @examples
#' library(igraph)
#' library(ggraph)

#' g <- sample_gnp(10,0.4)
#' coords <- layout_with_focus(g,v = 1)
#' coords
#' @export

layout_with_focus <- function(g,v,weights = NA,iter = 500,tol = 0.0001){
  if(!igraph::is.igraph(g)){
    stop("g must be an igraph object")
  }
  if(missing(v)){
    stop('argument "v" is missing with no default.')
  }
  comps <- igraph::components(g,"weak")
  if(comps$no>1){
    stop("g must be a connected graph.")
  }
  if (exists(".Random.seed", .GlobalEnv)){
    oldseed <- .GlobalEnv$.Random.seed
  }
  else{
    oldseed <- NULL
  }
  set.seed(42)  #stress is deterministic and produces same result up to translation. This keeps the layout fixed
  on.exit(restore_seed(oldseed))

  n <- igraph::vcount(g)
  D <- igraph::distances(g,weights = weights)
  W <- 1/D^2
  diag(W) <- 0

  Z <- matrix(0,n,n)
  Z[v,] <- Z[,v] <- 1
  Z <- W*Z


  rmat <- matrix(stats::runif(n*2,-0.1,0.1),n,2)
  xinit <- igraph::layout_with_mds(g) + rmat

  tseq <- seq(0,1,0.1)
  x <- stress_focus(xinit,W,D,Z,tseq,iter,tol)

  offset <- x[v,]
  x <- t(apply(x,1,function(x) x-offset))
  return(list(xy=x,distance=D[,v]))
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' radial centrality layout
#'
#' @description arranges nodes in concentric circles according to a centrality index.
#'
#' @name layout_centrality
#' @param g igraph object
#' @param cent centrality scores
#' @param scale logical. should centrality scores be scaled to \eqn{[0,100]}? (Default: TRUE)
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @param tseq numeric vector. increasing sequence of coefficients to combine regular stress and constraint stress. See details.
#' @details The function optimizes a convex combination of regular stress and a constrained stress function which forces
#' nodes to be arranged on concentric circles. The vector `tseq` is the sequence of parameters used for the convex combination.
#' In iteration i of the algorithm \eqn{tseq[i]} is used to combine regular and constraint stress as \eqn{(1-tseq[i])*stress_{regular}+tseq[i]*stress_{constraint}}. The sequence must be increasing, start at zero and end at one. The default setting should be a good choice for most graphs.
#'
#' The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' 'ggraph' natively supports the layout.
#' @return matrix of xy coordinates
#' @references Brandes, U., & Pich, C. (2011). More flexible radial layout. Journal of Graph Algorithms and Applications, 15(1), 157-173.
#' @examples
#' library(igraph)
#' library(ggraph)
#'
#' g <- sample_gnp(10,0.4)
#' \dontrun{
#' ggraph(g,layout="centrality",centrality = closeness(g))+
#'   draw_circle(use = "cent")+
#'   geom_edge_link0()+
#'   geom_node_point(shape = 21,fill = "grey25",size = 5)+
#'   theme_graph()+
#'   coord_fixed()
#'}

#' @export
#'
layout_with_centrality <- function(g,cent,scale = TRUE,iter = 500,tol = 0.0001,tseq = seq(0,1,0.2)){
  if(!igraph::is.igraph(g)){
    stop("g must be an igraph object")
  }
  comps <- igraph::components(g,"weak")
  if(comps$no>1){
    stop("g must be connected")
  }
  if (exists(".Random.seed", .GlobalEnv)){
    oldseed <- .GlobalEnv$.Random.seed
  }
  else{
    oldseed <- NULL
  }
  set.seed(42)  #stress is deterministic and produces same result up to translation. This keeps the layout fixed
  on.exit(restore_seed(oldseed))

  n <- igraph::vcount(g)
  if(scale){
    cent <- scale_to_100(cent)
  }
  r <- unname(igraph::diameter(g)/2 * (1 - ((cent-min(cent))/(max(cent)-min(cent)+1))))

  D <- igraph::distances(g,weights = NA)
  W <- 1/D^2
  diag(W) <- 0

  rmat <- matrix(stats::runif(n*2,-0.1,0.1),n,2)
  xinit <- igraph::layout_with_mds(g) + rmat

  x <- stress_major(xinit,W,D,iter,tol)
  x <- stress_radii(x,W,D,r,tseq)

  # move highest cent to 0,0
  idx <- which.max(cent)[1]
  offset <- x[idx,]

  x <- t(apply(x,1,function(x) x-offset))
  if(scale){
    radii_new <- round(100-cent,8)
    angles <- apply(x,1,function(y) atan2(y[2],y[1]))
    x <- cbind(radii_new *cos(angles),radii_new*sin(angles))
  } else{
    radii_new <- round(max(cent)-cent,8)
    angles <- apply(x,1,function(y) atan2(y[2],y[1]))
    x <- cbind(radii_new*cos(angles),radii_new*sin(angles))
    }
  x
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#' constrained stress layout
#'
#' @name layout_constrained_stress
#' @description force-directed graph layout based on stress majorization with variable constrained
#' @param g igraph object
#' @param coord numeric vector. fixed coordinates for dimension specified in `fixdim`.
#' @param fixdim string. which dimension should be fixed. Either "x" or "y".
#' @param weights possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @param mds should an MDS layout be used as initial layout (default: TRUE)
#' @param bbox constrain dimension of output. Only relevant to determine the placement of disconnected graphs
#' @details Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight).
#'
#' The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' 'ggraph' natively supports the layout.
#' @seealso [layout_constrained_stress3D]
#' @return matrix of xy coordinates
#' @references Gansner, E. R., Koren, Y., & North, S. (2004). Graph drawing by stress majorization. *In International Symposium on Graph Drawing* (pp. 239-250). Springer, Berlin, Heidelberg.
#' @export
layout_with_constrained_stress <- function(g,coord,fixdim="x",weights = NA,
                                           iter = 500,tol = 0.0001,mds = TRUE,bbox = 30){
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  if (exists(".Random.seed", .GlobalEnv)){
    oldseed <- .GlobalEnv$.Random.seed
  }
  else{
    oldseed <- NULL
  }
  set.seed(42)  #stress is deterministic and produces same result up to translation. This keeps the layout fixed
  on.exit(restore_seed(oldseed))

  fixdim <- match.arg(fixdim,c("x","y"))
  fixdim <- ifelse(fixdim=="x",1,2)

  if(missing(coord)){
    stop('"coord" is missing with no default.')
  }
  comps <- igraph::components(g,"weak")
  if(comps$no>1){
    stop('g must be connected')
  }

  if(igraph::vcount(g)==1){
    x <- matrix(c(0,0),1,2)
  } else{
    D <- igraph::distances(g,weights = weights)
    W <- 1/D^2
    diag(W) <- 0
    n <- igraph::vcount(g)
    if(!mds){
      xinit <- matrix(stats::runif(n*2,0,1),n,2)
      xinit[,fixdim] <- coord
    } else{
      rmat <- matrix(stats::runif(n*2,-0.1,0.1),n,2)
      if(igraph::vcount(g)<=100){
        xinit <- igraph::layout_with_mds(g) + rmat
      } else{
        xinit <- layout_with_pmds(g,D = D[,sample(1:(igraph::vcount(g)),100)]) + rmat
      }
      xinit[,fixdim] <- coord
    }
    x <- constrained_stress_major(xinit,fixdim,W,D,iter,tol)
  }
  x
}

#' constrained stress layout in 3D
#'
#' @name layout_constrained_stress3D
#' @description force-directed graph layout based on stress majorization with variable constrained in 3D
#' @param g igraph object
#' @param coord numeric vector. fixed coordinates for dimension specified in `fixdim`.
#' @param fixdim string. which dimension should be fixed. Either "x", "y" or "z".
#' @param weights possibly a numeric vector with edge weights. If this is NULL and the graph has a weight edge attribute, then the attribute is used. If this is NA then no weights are used (even if the graph has a weight attribute). By default, weights are ignored. See details for more.
#' @param iter number of iterations during stress optimization
#' @param tol stopping criterion for stress optimization
#' @param mds should an MDS layout be used as initial layout (default: TRUE)
#' @param bbox constrain dimension of output. Only relevant to determine the placement of disconnected graphs
#' @details Be careful when using weights. In most cases, the inverse of the edge weights should be used to ensure that the endpoints of an edges with higher weights are closer together (weights=1/E(g)$weight).
#'
#' This function does not come with direct support for igraph or ggraph.
#'
#' @seealso [layout_constrained_stress]
#' @return matrix of xyz coordinates
#' @references Gansner, E. R., Koren, Y., & North, S. (2004). Graph drawing by stress majorization. *In International Symposium on Graph Drawing* (pp. 239-250). Springer, Berlin, Heidelberg.
#' @export
layout_with_constrained_stress3D <- function(g,coord,fixdim="x",weights = NA,
                                             iter = 500,tol = 0.0001,mds = TRUE,bbox = 30){
  if (!igraph::is_igraph(g)) {
    stop("Not a graph object")
  }
  if (exists(".Random.seed", .GlobalEnv)){
    oldseed <- .GlobalEnv$.Random.seed
  }
  else{
    oldseed <- NULL
  }
  set.seed(42)
  on.exit(restore_seed(oldseed))
  fixdim <- match.arg(fixdim,c("x","y","z"))
  fixdim <- ifelse(fixdim=="x",1,ifelse(fixdim=="y",2,3))

  if(missing(coord)){
    stop('"coord" is missing with no default.')
  }
  comps <- igraph::components(g,"weak")
  if(comps$no>1){
    stop('g must be connected')
  }

  if(igraph::vcount(g)==1){
    x <- matrix(c(0,0,0),1,3)
  } else{
    D <- igraph::distances(g,weights = weights)
    W <- 1/D^2
    diag(W) <- 0
    n <- igraph::vcount(g)
    if(!mds){
      xinit <- matrix(stats::runif(n*3,0,1),n,3)
      xinit[,fixdim] <- coord
    } else{
      n <- igraph::vcount(g)
      pivs <- sample(1:n,min(c(50,n)))
      D1 <- D[,pivs]
      cmean <- colMeans(D1^2)
      rmean <- rowMeans(D1^2)
      Dmat <- D1^2-outer(rmean,cmean, function(x,y) x+y)+mean(D1^2)
      sl2 <- svd(Dmat)
      rmat <- matrix(stats::runif(n*3,-0.1,0.1),n,3)
      xinit <- (Dmat%*%sl2$v[,1:3]) + rmat
      row.names(xinit) <- NULL

      xinit[,fixdim] <- coord
    }
    x <- constrained_stress_major3D(xinit,fixdim,W,D,iter,tol)
  }
  x
}

#-------------------------------------------------------------------------------#
# helper functions ----
#-------------------------------------------------------------------------------#

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

scale_to_100 <- function(x){
  a <- min(x)
  b <- max(x)
  100/(b-a) * x -100/(b-a)*a
}

interpolate_cent <- function(cent,x){
  a <- min(cent)
  b <- max(cent)
  alpha <- 100/(b-a)
  beta <- -100/(b-a)*a
  (x-beta)/alpha
}

restore_seed <- function(oldseed){
  if (!is.null(oldseed))
    .GlobalEnv$.Random.seed <- oldseed
  else
    rm(".Random.seed", envir = .GlobalEnv)
}

#' @useDynLib graphlayouts, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
