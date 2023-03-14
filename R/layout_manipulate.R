#' @title manipulate layout
#' @description functions to manipulate an existing layout
#'
#' @param xy graph layout
#' @param angle angle for rotation
#' @param axis mirror horizontal or vertical
#' @name layout_manipulate
#' @details These functions are mostly useful for deterministic layouts such as [layout_with_stress]
#' @return manipulated matrix of xy coordinates
#' @examples
#' library(igraph)
#' g <- sample_gnp(50,0.3)
#'
#' xy <- layout_with_stress(g)
#'
#' #rotate 90 degrees
#' xy <- layout_rotate(xy,90)
#'
#' # flip horizontally
#' xy <- layout_mirror(xy,"horizontal")
#'
#' @author David Schoch
NULL

#' @rdname layout_manipulate
#' @export

layout_rotate <- function(xy,angle){
  if(!"matrix"%in%class(xy)){
    stop("xy must be a matrix")
  }
  if(!is.numeric(angle)){
    stop("angle must be numeric")
  }
  radians <- angle * pi / 180
  s <- sin(radians)
  c <- cos(radians)

  cbind(xy[ ,1] * c - xy[ ,2] * s,xy[ ,1] * s + xy[ ,2] * c)

}

#' @rdname layout_manipulate
#' @export

layout_mirror <- function(xy,axis = "vertical"){
  if(!"matrix"%in%class(xy)){
    stop("xy must be a matrix")
  }
  axis <- match.arg(axis,c("horizontal","vertical"))
  if(axis=="horizontal"){
    middle <- mean(xy[,2])
    ynew <- middle - (xy[,2] - middle)
    xnew <- xy[,1]
  } else if(axis=="vertical"){
    middle <- mean(xy[,1])
    xnew <- middle - (xy[,1] - middle)
    ynew <- xy[,2]
  }
  xynew <- cbind(xnew,ynew)
  colnames(xynew) <- NULL
  xynew
}
