# stress majorization layout in 3D

force-directed graph layout based on stress majorization in 3D.

## Usage

``` r
layout_with_stress3D(
  g,
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

## Value

matrix of xyz coordinates

## Details

Be careful when using weights. In most cases, the inverse of the edge
weights should be used to ensure that the endpoints of an edges with
higher weights are closer together (weights=1/E(g)\$weight).

## References

Gansner, E. R., Koren, Y., & North, S. (2004). Graph drawing by stress
majorization. *In International Symposium on Graph Drawing* (pp.
239-250). Springer, Berlin, Heidelberg.

## See also

[layout_stress](https://schochastics.github.io/graphlayouts/reference/layout_stress.md)
