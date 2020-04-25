#' backbone graph layout
#' @description emphasizes a hidden group structure if it exists in the graph. Calculates a layout for a sparsified network only including the most embedded edges. Deleted edges are added back after the layout is calculated.
#' @name layout_backbone
#' @param g igraph object
#' @param keep fraction of edges to keep during backbone calculation
#' @param backbone logical. Return edge ids of the backbone (Default: TRUE)
#' @details
#' The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
#' 'ggraph' natively supports the layout.
#' @return list of xy coordinates and vector of edge ids included in the backbone
#' @examples
#'  library(igraph)
#'
#'  g <- sample_islands(9,20,0.4,9)
#'  g <- simplify(g)
#'  V(g)$grp <- as.character(rep(1:9,each=20))
#'  bb <- layout_as_backbone(g,keep=0.4)
#'
#'  # add backbone links as edge attribute
#'  E(g)$col <- FALSE
#'  E(g)$col[bb$backbone] <- TRUE
#'
#'
#' @references Nocaj, A., Ortmann, M., & Brandes, U. (2015). Untangling the hairballs of multi-centered, small-world online social media networks. Journal of Graph Algorithms and Applications: JGAA, 19(2), 595-618.
#' @export
#'

layout_as_backbone <- function(g,keep=0.2,backbone = TRUE){
  if(!requireNamespace("oaqc", quietly = TRUE)){
    stop("oaqc is needed for this function to work. Please install it.", call. = FALSE)
  }
  if(igraph::any_multiple(g)){
    stop("backbone layout does not work with multiple edges.")
  }
  if(igraph::is_directed(g)){
    stop("backbone layout does not work with directed edges.")
  }
  if(any(igraph::is.loop(g))){
    stop("backbone layout does not work with loops.")
  }
  #weighting ----
  orbs <- oaqc::oaqc(igraph::get.edgelist(g,names = FALSE)-1, non_ind_freq = T)
  e11 <- orbs$e_orbits_non_ind[ ,11]

  qu <- rep(0, igraph::vcount(g))
  el <- igraph::get.edgelist(g, names = FALSE)
  el <- cbind(el, e11)
  for(e in 1:nrow(el)){
    qu[el[e,1]] <- qu[el[e,1]] + el[e,3]
    qu[el[e,2]] <- qu[el[e,2]] + el[e,3]
  }
  w <- apply(el,1,function(x) x[3]/sqrt(qu[x[1]] * qu[x[2]]))

  w[is.na(w)] <- 0
  w[is.infinite(w)] <- 0
  igraph::E(g)$weight <- w

  #reweighting -----
  w <- max_prexif_jaccard(g)
  igraph::E(g)$weight <- w

  # umst ----
  g_umst <- umst(g)

  #filtering ----
  igraph::E(g)$bone <- w>=sort(w,decreasing=T)[ceiling(igraph::ecount(g) * keep)]
  g_bone <- igraph::graph_from_edgelist(el[igraph::E(g)$bone,1:2],directed = F)
  g_lay <- igraph::simplify(igraph::graph.union(g_umst,g_bone))
  if(backbone){
    bb <- backbone_edges(g,g_lay)
  } else {
    bb <- NULL
  }
  xy <- layout_with_stress(g_lay)
  list(xy = xy,backbone = bb)
}

#-------------------------------------------------------------------------------
# helper functions
#-------------------------------------------------------------------------------

umst <- function(g){
  el <- igraph::get.edgelist(g,names = FALSE)
  el <- cbind(el,igraph::E(g)$weight)
  el <- el[order(el[,3], decreasing = TRUE),]
  el <- cbind(el,rank(-el[,3]))
  vfind <- 1:igraph::vcount(g)
  el_un <- matrix(0,0,2)
  for(i in unique(el[,4])){
    el_tmp <- matrix(0,0,2)
    Bi <- which(el[,4]==i)
    for(e in Bi){
      u <- el[e,1]
      v <- el[e,2]
      if(vfind[u]!=vfind[v]){
        el_tmp <- rbind(el_tmp,c(u,v))
      }
    }
    if(nrow(el_tmp)==0){
      next()
    }
    for(eb in 1:nrow(el_tmp)){
      u <- el_tmp[eb,1]
      v <- el_tmp[eb,2]
      partu <- vfind[u]
      partv <- vfind[v]
      vfind[v] <- partu
      if(any(vfind==partv)){
        vfind[vfind==partv] <- partu
      }
    }
    el_un <- rbind(el_un,el_tmp)
  }
  return(igraph::simplify(igraph::graph_from_edgelist(el_un,directed=FALSE)))
}


backbone_edges <- function(g,g_lay){
  tmp <- rbind(igraph::get.edgelist(g_lay),igraph::get.edgelist(g,names = FALSE))
  which(duplicated(tmp))-igraph::ecount(g_lay)
}

max_prexif_jaccard <- function(g){

  if("name"%in%igraph::vertex_attr_names(g)){
    g <- igraph::delete_vertex_attr(g,"name")
  }
  el_tbl <- igraph::as_data_frame(g,"edges")

  N_ranks <- lapply(1:igraph::vcount(g),get_rank,el_tbl=el_tbl)
  el <- igraph::get.edgelist(g,names = F)
  new_w <- reweighting(el-1,N_ranks)
  new_w
}

get_rank <- function(el_tbl,u){
  Nu_idx <- el_tbl[["from"]]==u | el_tbl[["to"]]==u
  omega <- el_tbl[Nu_idx,"weight"]
  Nu <- setdiff(c(el_tbl[Nu_idx,"from"],el_tbl[Nu_idx,"to"]),u)
  r <- rank(-omega)
  r <- match(r, sort(unique(r)))-1
  Nru <- cbind(Nu-1,r)
  Nru[order(Nru[,2]),,drop = FALSE]
}
