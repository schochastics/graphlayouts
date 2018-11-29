#' Backbone layout graph layout
#'
#' @param g igraph object
#' @param keep fraction of edges to keep
#' @param backbone logical. Return edge ids of the backbone
#'
#' @return coordinates to be used layouting a graph
#' @export
#'

backbone_layout <- function(g,keep=0.2,backbone=T){
  orbs <- oaqc::oaqc(igraph::get.edgelist(g,names = F)-1)
  qu <- rowSums(orbs$n_orbits_ind[,17:20])
  el <- igraph::get.edgelist(g)
  el <- cbind(el,rowSums(orbs$e_orbits_ind[,11:14]))

  w <- apply(el,1,function(x) x[3]/sqrt(qu[x[1]]*qu[x[2]]))
  w[is.na(w)] <- 0
  w[is.infinite(w)] <- 0
  igraph::E(g)$weight <- w
  igraph::E(g)$bone <- w>=sort(w,decreasing=T)[ceiling(igraph::ecount(g)*keep)]
  g_umst <- umst(g)
  g_bone <- igraph::graph_from_edgelist(el[igraph::E(g)$bone,1:2],directed = F)
  g_lay <- igraph::simplify(igraph::graph.union(g_umst,g_bone))
  if(backbone){
    bb <- backbone_edges(g,g_lay)
  } else {
    bb <- NULL
  }
  xy <- smglr::stress_majorization(g_lay)
  list(xy=xy,backbone=bb)
}

#-------------------------------------------------------------------------------
# helper functions
#-------------------------------------------------------------------------------

umst <- function(g){
  el <- igraph::get.edgelist(g)
  el <- cbind(el,igraph::E(g)$weight)
  el <- el[order(el[,3],decreasing=T),]
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
    for(e in nrow(el_tmp)){
      u <- el_tmp[e,1]
      v <- el_tmp[e,2]
      partu <- vfind[u]
      partv <- vfind[v]
      vfind[v] <- partu
      vfind[vfind==partv] <- partu
    }
    el_un <- rbind(el_un,el_tmp)
  }
  return(igraph::simplify(igraph::graph_from_edgelist(el_un,directed=F)))
}


backbone_edges <- function(g,g_lay){
  tmp <- rbind(igraph::get.edgelist(g_lay),igraph::get.edgelist(g))
  which(duplicated(tmp))-igraph::ecount(g_lay)
}

