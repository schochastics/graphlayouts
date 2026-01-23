# pivot MDS graph layout

similar to
[layout_with_mds](https://r.igraph.org/reference/layout_with_mds.html)
but uses only a small set of pivots for MDS. Considerably faster than
MDS and thus applicable for larger graphs.

## Usage

``` r
layout_with_pmds(g, pivots, weights = NA, D = NULL, dim = 2)

layout_igraph_pmds(g, pivots, weights = NA, D = NULL, circular)
```

## Arguments

- g:

  igraph object

- pivots:

  number of pivots

- weights:

  possibly a numeric vector with edge weights. If this is NULL and the
  graph has a weight edge attribute, then the attribute is used. If this
  is NA then no weights are used (even if the graph has a weight
  attribute). By default, weights are ignored. See details for more.

- D:

  precomputed distances from pivots to all nodes (if available, default:
  NULL)

- dim:

  dimensionality of layout (defaults to 2)

- circular:

  not used

## Value

matrix of coordinates

## Details

Be careful when using weights. In most cases, the inverse of the edge
weights should be used to ensure that the endpoints of an edges with
higher weights are closer together (weights=1/E(g)\$weight)

The layout_igraph\_\* function should not be used directly. It is only
used as an argument for plotting with 'igraph'. 'ggraph' natively
supports the layout.

## References

Brandes, U. and Pich, C. (2006). Eigensolver Methods for Progressive
Multidimensional Scaling of Large Data. In *International Symposium on
Graph Drawing* (pp. 42-53). Springer

## Author

David Schoch

## Examples

``` r
if (FALSE) { # \dontrun{
library(igraph)
library(ggraph)

g <- sample_gnp(1000, 0.01)

xy <- layout_with_pmds(g, pivots = 100)
} # }
```
