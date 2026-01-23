# sparse stress graph layout

stress majorization for larger graphs based on a set of pivot nodes.

## Usage

``` r
layout_with_sparse_stress(g, pivots, weights = NA, iter = 500)

layout_igraph_sparse_stress(g, pivots, weights = NA, iter = 500, circular)
```

## Arguments

- g:

  igraph object

- pivots:

  number of pivots

- weights:

  ignored

- iter:

  number of iterations during stress optimization

- circular:

  not used

## Value

matrix of xy coordinates

## Details

The layout_igraph\_\* function should not be used directly. It is only
used as an argument for plotting with 'igraph'. 'ggraph' natively
supports the layout.

## References

Ortmann, M. and Klimenta, M. and Brandes, U. (2016). A Sparse Stress
Model. https://arxiv.org/pdf/1608.08909.pdf

## Author

David Schoch

## Examples

``` r
if (FALSE) { # \dontrun{
library(igraph)
library(ggraph)

g <- sample_gnp(1000, 0.005)

ggraph(g, layout = "sparse_stress", pivots = 100) +
    geom_edge_link0(edge_colour = "grey66") +
    geom_node_point(shape = 21, fill = "grey25", size = 5) +
    theme_graph()
} # }
```
