# radial centrality layout

arranges nodes in concentric circles according to a centrality index.

## Usage

``` r
layout_with_centrality(
  g,
  cent,
  scale = TRUE,
  iter = 500,
  tol = 1e-04,
  tseq = seq(0, 1, 0.2)
)

layout_igraph_centrality(
  g,
  cent,
  scale = TRUE,
  iter = 500,
  tol = 1e-04,
  tseq = seq(0, 1, 0.2),
  circular
)
```

## Arguments

- g:

  igraph object

- cent:

  centrality scores

- scale:

  logical. should centrality scores be scaled to \\\[0,100\]\\?
  (Default: TRUE)

- iter:

  number of iterations during stress optimization

- tol:

  stopping criterion for stress optimization

- tseq:

  numeric vector. increasing sequence of coefficients to combine regular
  stress and constraint stress. See details.

- circular:

  not used

## Value

matrix of xy coordinates

## Details

The function optimizes a convex combination of regular stress and a
constrained stress function which forces nodes to be arranged on
concentric circles. The vector `tseq` is the sequence of parameters used
for the convex combination. In iteration i of the algorithm
\\tseq\[i\]\\ is used to combine regular and constraint stress as
\\(1-tseq\[i\])\*stress\_{regular}+tseq\[i\]\*stress\_{constraint}\\.
The sequence must be increasing, start at zero and end at one. The
default setting should be a good choice for most graphs.

The layout_igraph\_\* function should not be used directly. It is only
used as an argument for plotting with 'igraph'. 'ggraph' natively
supports the layout.

## References

Brandes, U., & Pich, C. (2011). More flexible radial layout. Journal of
Graph Algorithms and Applications, 15(1), 157-173.

## See also

[layout_centrality_group](https://schochastics.github.io/graphlayouts/reference/layout_centrality_group.md)

## Examples

``` r
library(igraph)

g <- sample_gnp(10, 0.4)
if (FALSE) { # \dontrun{
library(ggraph)
ggraph(g, layout = "centrality", centrality = closeness(g)) +
    draw_circle(use = "cent") +
    geom_edge_link0() +
    geom_node_point(shape = 21, fill = "grey25", size = 5) +
    theme_graph() +
    coord_fixed()
} # }
```
