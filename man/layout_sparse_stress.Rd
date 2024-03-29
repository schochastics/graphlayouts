% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/layout_large_graphs.R, R/layouts.R
\name{layout_sparse_stress}
\alias{layout_sparse_stress}
\alias{layout_with_sparse_stress}
\alias{layout_igraph_sparse_stress}
\title{sparse stress graph layout}
\usage{
layout_with_sparse_stress(g, pivots, weights = NA, iter = 500)

layout_igraph_sparse_stress(g, pivots, weights = NA, iter = 500, circular)
}
\arguments{
\item{g}{igraph object}

\item{pivots}{number of pivots}

\item{weights}{ignored}

\item{iter}{number of iterations during stress optimization}

\item{circular}{not used}
}
\value{
matrix of xy coordinates
}
\description{
stress majorization for larger graphs based on a set of pivot nodes.
}
\details{
The layout_igraph_* function should not be used directly. It is only used as an argument for plotting with 'igraph'.
'ggraph' natively supports the layout.
}
\examples{
\dontrun{
library(igraph)
library(ggraph)

g <- sample_gnp(1000, 0.005)

ggraph(g, layout = "sparse_stress", pivots = 100) +
    geom_edge_link0(edge_colour = "grey66") +
    geom_node_point(shape = 21, fill = "grey25", size = 5) +
    theme_graph()
}
}
\references{
Ortmann, M. and Klimenta, M. and Brandes, U. (2016). A Sparse Stress Model. https://arxiv.org/pdf/1608.08909.pdf
}
\author{
David Schoch
}
