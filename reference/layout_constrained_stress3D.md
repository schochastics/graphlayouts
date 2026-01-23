# constrained stress layout in 3D

force-directed graph layout based on stress majorization with variable
constrained in 3D

## Usage

``` r
layout_with_constrained_stress3D(
  g,
  coord,
  fixdim = "x",
  weights = NA,
  iter = 500,
  tol = 1e-04,
  mds = TRUE,
  bbox = 30
)
```

## Arguments

- g:

  igraph object

- coord:

  numeric vector. fixed coordinates for dimension specified in `fixdim`.

- fixdim:

  string. which dimension should be fixed. Either "x", "y" or "z".

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

## Value

matrix of xyz coordinates

## Details

Be careful when using weights. In most cases, the inverse of the edge
weights should be used to ensure that the endpoints of an edges with
higher weights are closer together (weights=1/E(g)\$weight).

This function does not come with direct support for igraph or ggraph.

## References

Gansner, E. R., Koren, Y., & North, S. (2004). Graph drawing by stress
majorization. *In International Symposium on Graph Drawing* (pp.
239-250). Springer, Berlin, Heidelberg.

## See also

[layout_constrained_stress](https://schochastics.github.io/graphlayouts/reference/layout_constrained_stress.md)
