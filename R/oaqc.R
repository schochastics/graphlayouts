#' This files is adapted from the archived R package oaqc
#' Coerce graph input.
#'
#' @param graph A matrix, data.frame or graph object.
#' @return Edge list matrix.
as.edge_list <- function(graph) {
    if (is.matrix(graph)) {
        if (length(dim(graph)) != 2 && ncol(graph) != 2) {
            stop(paste(
                "not coercable to edge list",
                "only matrices with 2 columns supported"
            ))
        }
        graph
    } else if (is.data.frame(graph)) {
        if (ncol(graph) != 2) {
            stop(paste(
                "not coercable to edge list",
                "only data frames with 2 columns supported"
            ))
        }
        data.matrix(graph)
    } else if (inherits(graph, "igraph") && requireNamespace("igraph", quietly = TRUE)) {
        igraph::as_edgelist(graph) - 1
    } else {
        stop(paste(
            "unrecognized graph type",
            "use matrix/data frame with 2 columns or igraph objects"
        ))
    }
}

#' Annotates the igraph object with orbit labels.
#'
#' @param graph Unmodified input graph.
#' @param orbits List with n_orbits, e_orbits matrices.
#' @param non_ind_freq A flag indicating whether non-induced frequencies have to be written or not.
#' @return \code{orbits} if the input is not an igraph, the annotated igraph
#' instead.
annotate_result <- function(graph, orbits, non_ind_freq) {
    if (inherits(graph, "igraph") && requireNamespace("igraph", quietly = T)) {
        for (i in seq_len(ncol(orbits$n_orbits_ind))) {
            graph <- igraph::set_vertex_attr(graph,
                paste("orbit_ind", i, sep = "_"),
                value = orbits$n_orbits_ind[, i]
            )
        }
        for (i in seq_len(col(orbits$e_orbits_ind))) {
            graph <- igraph::set_edge_attr(graph,
                paste("orbit_ind", i, sep = "_"),
                value = orbits$e_orbits_ind[, i]
            )
        }
        if (non_ind_freq) {
            for (i in seq_len(ncol(orbits$`n_orbits_non_ind`))) {
                graph <- igraph::set_vertex_attr(graph,
                    paste("orbit_non_ind", i, sep = "_"),
                    value = orbits$`n_orbits_non_ind`[, i]
                )
            }
            for (i in seq_len(ncol(orbits$`e_orbits_non_ind`))) {
                graph <- igraph::set_edge_attr(graph,
                    paste("orbit_non_ind", i, sep = "_"),
                    value = orbits$`e_orbits_non_ind`[, i]
                )
            }
        }
        return(graph)
    }
    return(orbits)
}

#' Calculates the orbit-aware quad census on an edge and node level, see
#' \code{vignette('oaqc')}.
#'
#' @param graph A matrix, data.frame or graph object.
#' @param non_ind_freq A flag indicating whether non-induced frequencies have to be returned or not.
#' @param file Name (and location) of the file to be written.
#' @return orbit-aware quad census on a node and edge level. Consult
#' \code{vignette('oaqc')} to see the correspondence between orbit and quad.
#' @examples
#' k4 <- data.frame(
#'     source = c(0, 0, 0, 1, 1, 2),
#'     target = c(1, 2, 3, 2, 3, 3)
#' )
#'
#' k4orbits <- oaqc(k4, non_ind_freq = TRUE)
#' print(k4orbits)
#' @export
oaqc <- function(graph, non_ind_freq = F, file = "") {
    edges <- as.edge_list(graph)
    orbits <- .Call(
        entry, as.integer(max(edges) + 1), as.integer(edges),
        as.logical(non_ind_freq), as.character(file)
    )
    return(annotate_result(graph, orbits, as.logical(non_ind_freq)))
}
