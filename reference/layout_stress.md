# stress majorization layout

force-directed graph layout based on stress majorization. Similar to
Kamada-Kawai, but generally faster and with better results.

## Usage

``` r
layout_with_stress(
  g,
  weights = NA,
  iter = 500,
  tol = 1e-04,
  mds = TRUE,
  bbox = 30
)

layout_igraph_stress(
  g,
  weights = NA,
  iter = 500,
  tol = 1e-04,
  mds = TRUE,
  bbox = 30,
  circular
)
```

## Arguments

- g:

  igraph object

- weights:

  possibly a numeric vector with edge weights. If this is NULL and the
  graph has a weight edge attribute, then the attribute is used. If this
  is NA then no weights are used (even if the graph has a weight
  attribute). By default, weights are ignored. See details for more.

- iter:

  number of iterations during stress optimization

- tol:

  stopping criterion for stress optimization

- mds:

  should an MDS layout be used as initial layout (default: TRUE)

- bbox:

  width of layout. Only relevant to determine the placement of
  disconnected graphs

- circular:

  not used

## Value

matrix of xy coordinates

## Details

Be careful when using weights. In most cases, the inverse of the edge
weights should be used to ensure that the endpoints of an edges with
higher weights are closer together (weights=1/E(g)\$weight).

The layout_igraph\_\* function should not be used directly. It is only
used as an argument for plotting with 'igraph'. 'ggraph' natively
supports the layout.

## References

Gansner, E. R., Koren, Y., & North, S. (2004). Graph drawing by stress
majorization. *In International Symposium on Graph Drawing* (pp.
239-250). Springer, Berlin, Heidelberg.

## See also

[layout_stress3D](https://schochastics.github.io/graphlayouts/reference/layout_stress3D.md)

## Examples

``` r
library(igraph)
set.seed(665)

g <- sample_pa(100, 1, 1, directed = FALSE)

# calculate layout manually
xy <- layout_with_stress(g)

# use it with ggraph
if (FALSE) { # \dontrun{
library(ggraph)
ggraph(g, layout = "stress") +
    geom_edge_link0(edge_width = 0.2, colour = "grey") +
    geom_node_point(col = "black", size = 0.3) +
    theme_graph()
} # }
```
