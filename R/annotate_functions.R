#' Draw concentric circles
#'
#' @param col color of circles
#' @param use one of 'focus' or 'cent'
#' @param max.circle if use = 'focus' specifies the number of circles to draw
#' @details this function is best used with a concentric layout such as [layout_with_focus] and [layout_with_centrality].
#' @return concentric circles around origin
#' @examples
#' library(igraph)
#' library(ggraph)

#' g <- sample_gnp(10,0.4)
#'
#' \dontrun{
#' ggraph(g,layout = "centrality",centrality = degree(g))+
#'   draw_circle(use = "cent")+
#'   geom_edge_link()+
#'   geom_node_point(shape = 21,fill = "grey25",size = 5)+
#'   theme_graph()+
#'   coord_fixed()
#'}
#' @export
#'

draw_circle <- function(col = "#00BFFF",use="focus",max.circle){
  if(!requireNamespace("ggplot2", quietly = TRUE)){
    stop("ggplot2 needed for this function to work. Please install it.", call. = FALSE)
  }
  if(!use%in%c("focus","cent")){
    stop("use must be one of 'focus' or 'cent'")
  }
  if(use=="focus" & missing(max.circle)){
    stop("max.circle missing. Should be set to the max distance from focal node.")
  }
  dat <- data.frame()
  if(use=="cent"){
    for(d in seq(0,100,20)*2){
      tmp <- as.data.frame(circleFun(center=c(0,0),diameter = d))
      tmp[["grp"]] <- d
      dat <- rbind(dat,tmp)
    }
  } else if(use=="focus"){
    for(d in 1:max.circle){
      tmp <- as.data.frame(circleFun(center=c(0,0),diameter = 2*d))
      tmp[["grp"]] <- d
      dat <- rbind(dat,tmp)
    }
  }
  circs <- ggplot2::geom_path(data=dat,ggplot2::aes_(x = ~x,y = ~y,group = ~grp),col=col,alpha=0.5)
  return(circs)
}


#' annotate concentric circles
#'
#' @param cent centrality scores used for layout
#' @param col color of text
#' @param format either empty string or 'scientific'
#' @param pos position of text ('top' or 'bottom')
#' @param text_size font size for annotations
#' @details this function is best used with [layout_with_centrality] together with [draw_circle].
#' @return annotated concentric circles around origin
#' @examples
#'library(igraph)
#'library(ggraph)
#'
#'g <- sample_gnp(10,0.4)
#'\dontrun{
#'ggraph(g,layout = "centrality",centrality = closeness(g))+
#'  draw_circle(use = "cent")+
#'  annotate_circle(closeness(g),pos = "bottom",format = "scientific")+
#'  geom_edge_link()+
#'  geom_node_point(shape=21,fill="grey25",size=5)+
#'  theme_graph()+
#'  coord_fixed()
#'}
#' @export
#'
annotate_circle <- function(cent,col = "#00BFFF",format="",pos="top",text_size=3){
  if(!requireNamespace("ggplot2", quietly = TRUE)){
    stop("ggplot2 needed for this function to work. Please install it.", call. = FALSE)
  }
  if(length(cent)==1){
    cent <- seq(1,cent,1)
    dat_annot <- data.frame(y=seq(0,100,20),x=0,val=interpolate_cent(cent,seq(0,100,20)))
    dat_annot[["val"]] <- round(dat_annot[["val"]],8)
  } else{
    dat_annot <- data.frame(y=100-seq(0,100,20),x=0,val=interpolate_cent(cent,seq(0,100,20)))
    dat_annot[["val"]] <- round(dat_annot[["val"]],8)
  }
  vju <- 0
  if(format=="scientific"){
    dat_annot[["val"]] <- format(dat_annot[["val"]],scientific = TRUE)
  }
  if(pos=="bottom"){
    dat_annot[["y"]] <- -dat_annot[["y"]]
    vju <- 1
  }

  circs <-ggplot2::geom_text(data=dat_annot,ggplot2::aes_(x = ~x,y = ~y,label = ~val),
                             col=col,vjust=vju,size=text_size)
  return(circs)
}


#helper -----
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

interpolate_cent <- function(cent,x){
  a <- min(cent)
  b <- max(cent)
  alpha <- 100/(b-a)
  beta <- -100/(b-a)*a
  (x-beta)/alpha
}
