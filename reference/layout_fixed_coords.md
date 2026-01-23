# Layout with fixed coordinates

force-directed graph layout based on stress majorization with fixed
coordinates for some nodes

## Usage

``` r
layout_with_fixed_coords(
  g,
  coords,
  weights = NA,
  iter = 500,
  tol = 1e-04,
  mds = TRUE,
  bbox = 30
)

layout_igraph_fixed_coords(
  g,
  coords,
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

- coords:

  numeric n x 2 matrix, where n is the number of nodes. values are
  either NA or fixed coordinates. coordinates are only calculated for
  the NA values.

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

  constrain dimension of output. Only relevant to determine the
  placement of disconnected graphs

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

## See also

[layout_constrained_stress](https://schochastics.github.io/graphlayouts/reference/layout_constrained_stress.md)

## Examples

``` r
library(igraph)
set.seed(12)
g <- sample_bipartite(10, 5, "gnp", 0.5)
#> Warning: `sample_bipartite()` was deprecated in igraph 2.2.0.
#> ℹ Please use `sample_bipartite_gnp()` instead.
fxy <- cbind(c(rep(0, 10), rep(1, 5)), NA)
xy <- layout_with_fixed_coords(g, fxy)
```
