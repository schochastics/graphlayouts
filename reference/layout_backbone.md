# backbone graph layout

emphasizes a hidden group structure if it exists in the graph.
Calculates a layout for a sparsified network only including the most
embedded edges. Deleted edges are added back after the layout is
calculated.

## Usage

``` r
layout_as_backbone(g, keep = 0.2, backbone = TRUE)

layout_igraph_backbone(g, keep = 0.2, backbone = TRUE, circular)
```

## Arguments

- g:

  igraph object

- keep:

  fraction of edges to keep during backbone calculation

- backbone:

  logical. Return edge ids of the backbone (Default: TRUE)

- circular:

  not used

## Value

list of xy coordinates and vector of edge ids included in the backbone

## Details

The layout_igraph\_\* function should not be used directly. It is only
used as an argument for plotting with 'igraph'. 'ggraph' natively
supports the layout.

## References

Nocaj, A., Ortmann, M., & Brandes, U. (2015). Untangling the hairballs
of multi-centered, small-world online social media networks. Journal of
Graph Algorithms and Applications: JGAA, 19(2), 595-618.

## Examples

``` r
library(igraph)
if (FALSE) { # \dontrun{
g <- sample_islands(9, 20, 0.4, 9)
g <- simplify(g)
V(g)$grp <- as.character(rep(1:9, each = 20))
bb <- layout_as_backbone(g, keep = 0.4)

# add backbone links as edge attribute
E(g)$col <- FALSE
E(g)$col[bb$backbone] <- TRUE
} # }
```
