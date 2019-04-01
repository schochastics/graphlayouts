#' Stress majorization graph layout
#'
#' @name stress_layout
#' @param g igraph object
#' @param iter number of iterations
#' @param tol stopping criterion
#' @param mds should an MDS layout be used as initial layout (default: TRUE)
#' @param bbox constrain dimension of output
#' @details the layout_igraph_* function should not be used directly. It is only used as an argument for 'ggraph'.
#' @return coordinates to be used layouting a graph
#' @references Gansner, E. R., Koren, Y., & North, S. (2004). Graph drawing by stress majorization. In International Symposium on Graph Drawing (pp. 239-250). Springer, Berlin, Heidelberg.
#' @examples
#' library(igraph)
#' library(ggraph)
#' set.seed(665)
#'
#' g <- sample_pa(100,1,1,directed = FALSE)
#'
#' #calculate layout manualy
#' xy <- layout_with_stress(g)
#' #use it with ggraph
#' ggraph(g,layout="stress")+
#'   geom_edge_link(width=0.2,colour="grey")+
#'   geom_node_point(col="black",size=0.3)+
#'   theme_graph()
#' @export
#'
layout_with_stress <- function(g,iter=500,tol=0.0001,mds=TRUE,bbox=50){
  if(!igraph::is.igraph(g)){
    stop("g must be an igraph object")
  }
  comps <- igraph::components(g,"weak")
  if(comps$no>1){
    lg <- list()
    node_order <- c()
    for (i in 1:comps$no){
      sg <- igraph::induced_subgraph(g,comps$membership==i)
      n <- igraph::vcount(sg)
      node_order <- c(node_order,which(comps$membership==i))
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
    x <- x[order(node_order),]

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
# focal layout
#-------------------------------------------------------------------------------
#' Focal layout
#'
#' @description puts a focal node in the center and arrange other nodes in concentric circles according to distances.
#'
#' @name focal_layout
#' @param g igraph object
#' @param v focal node to be placed in the center
#' @param iter number of iterations
#' @param tol stopping criterion
#' @details the layout_igraph_* function should not be used directly. It is only used as an argument for 'ggraph'.
#' @return coordinates to be used layouting a graph
#' @references Brandes, U., & Pich, C. (2011). More flexible radial layout. Journal of Graph Algorithms and Applications, 15(1), 157-173.
#' @examples
#' library(igraph)
#' library(ggraph)

#' g <- sample_gnp(10,0.4)

#' ggraph(g,layout = "focus",v = 1)+
#'   draw_circle(use = "focus", max.circle = max(distances(g,1)))+
#'   geom_edge_link()+
#'   geom_node_point(shape = 21,fill = "grey25",size = 5)+
#'   theme_graph()+
#'   coord_fixed()
#' @export
layout_with_focus <- function(g,v,iter=500,tol=0.0001){
  if(!igraph::is.igraph(g)){
    stop("g must be an igraph object")
  }
  if(missing(v)){
    stop("no focal node provided")
  }
  comps <- igraph::components(g,"weak")
  if(comps$no>1){
    stop("g must be connected")
  }
  n <- igraph::vcount(g)
  D <- igraph::distances(g,weights = NA)
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
  x
}

#-------------------------------------------------------------------------------
# centrality layout
#-------------------------------------------------------------------------------
#' Centrality layout
#'
#' @description arranges nodes in concentric circles according to a centrality index
#'
#' @name centrality_layout
#' @param g igraph object
#' @param cent centrality scores
#' @param scale scale centrality between 0 and 100?
#' @param iter number of iterations
#' @param tol stopping criterion
#' @param tseq transition steps
#' @details the layout_igraph_* function should not be used directly. It is only used as an argument for 'ggraph'.
#' @return coordinates to be used layouting a graph
#' @references Brandes, U., & Pich, C. (2011). More flexible radial layout. Journal of Graph Algorithms and Applications, 15(1), 157-173.
#' @examples
#' library(igraph)
#' library(ggraph)

#' g <- sample_gnp(10,0.4)

#' ggraph(g,layout="centrality",cent=closeness(g))+
#'   draw_circle(use = "cent")+
#'   geom_edge_link()+
#'   geom_node_point(shape=21,fill="grey25",size=5)+
#'   theme_graph()+
#'   coord_fixed()

#' @export
#'
layout_with_centrality <- function(g,cent,scale=T,iter=500,tol=0.0001,tseq=seq(0,1,0.2)){
  if(!igraph::is.igraph(g)){
    stop("g must be an igraph object")
  }
  comps <- igraph::components(g,"weak")
  if(comps$no>1){
    stop("g must be connected")
  }
  n <- igraph::vcount(g)
  if(scale){
    cent <- scale_to_100(cent)
  }
  r <- unname(igraph::diameter(g)/2 * (1 - ((cent-min(cent))/(max(cent)-min(cent)+1))))

  D <- igraph::distances(g,weights = NA)
  W <- 1/D^2
  diag(W) <- 0

  # D <- rbind(cbind(D,r),c(r,0))
  # W <- rbind(cbind(W,0),0)
  #
  # Z <- matrix(0,n+1,n+1)
  # Z[(n+1),1:n] <- Z[1:n,(n+1)] <- 1/(r^2)


  # xinit <- matrix(stats::runif((n+1)*2,0,1),n+1,2)
  rmat <- matrix(stats::runif(n*2,-0.1,0.1),n,2)
  xinit <- igraph::layout_with_mds(g) + rmat
  # xinit <- rbind(xinit,c(0,0))
  x <- stress_major(xinit,W,D,iter,tol)
  # x[n+1,] <- apply(x,2,mean)
  # tseq <- seq(0,1,0.1)
  x <- stress_radii(x,W,D,r,tseq)
  # x <- stress_focus(xinit,W,D,Z,tseq,iter,tol)
  # x <- x[1:n,]

  #move highest cent to 0,0
  idx <- which.max(cent)[1]
  offset <- x[idx,]
  # x[idx,] <- c(0,0)
  x <- t(apply(x,1,function(x) x-offset))
  if(scale){
    radii_new <- round(100-cent,8)
    angles <- apply(x,1,function(y) atan2(y[2],y[1]))
    x <- cbind(radii_new*cos(angles),radii_new*sin(angles))
  } else{
    radii_new <- round(max(cent)-cent,8)
    angles <- apply(x,1,function(y) atan2(y[2],y[1]))
    x <- cbind(radii_new*cos(angles),radii_new*sin(angles))
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

#' @useDynLib graphlayouts
#' @importFrom Rcpp sourceCpp
NULL
