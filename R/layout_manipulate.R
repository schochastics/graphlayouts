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
  radians <- angle * pi / 180
  s <- sin(radians)
  c <- cos(radians)

  xnew <- xy[ ,1] * c - xy[ ,2] * s
  ynew <- xy[ ,1] * s + xy[ ,2] * c

  cbind(xnew,ynew)

}

#' @rdname layout_manipulate
#' @export

layout_mirror <- function(xy,axis = "vertical"){
  if(grepl(axis,"vertical")){
    axis <- "vertical"
  } else if(grepl(axis,"horizontal")){
    axis <- "horizontal"
  } else{
    stop("axis must be vertical, horizontal or an abbreviation of the two.")
  }
  if(axis=="horizontal"){
    middle <- mean(xy[,2])
    ynew <- middle - (xy[,2] - middle)
    xnew <- xy[,1]
  } else if(axis=="vertical"){
    middle <- mean(xy[,1])
    xnew <- middle - (xy[,1] - middle)
    ynew <- xy[,2]
  }
  cbind(xnew,ynew)
}
