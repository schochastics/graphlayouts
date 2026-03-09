# radial focus layout

arrange nodes in concentric circles around a focal node according to
their distance from the focus.

## Usage

``` r
layout_with_focus(g, v, weights = NA, iter = 500, tol = 1e-04)

layout_igraph_focus(g, v, weights = NA, iter = 500, tol = 1e-04, circular)
```

## Arguments

- g:

  igraph object

- v:

  id of focal node to be placed in the center

- weights:

  possibly a numeric vector with edge weights. If this is NULL and the
  graph has a weight edge attribute, then the attribute is used. If this
  is NA then no weights are used (even if the graph has a weight
  attribute). By default, weights are ignored. See details for more.

- iter:

  number of iterations during stress optimization

- tol:

  stopping criterion for stress optimization

- circular:

  not used

## Value

a list containing xy coordinates and the distances to the focal node

## Details

Be careful when using weights. In most cases, the inverse of the edge
weights should be used to ensure that the endpoints of an edges with
higher weights are closer together (weights=1/E(g)\$weight).

## References

Brandes, U., & Pich, C. (2011). More flexible radial layout. *Journal of
Graph Algorithms and Applications*, 15(1), 157-173.

## See also

[layout_focus_group](https://schochastics.github.io/graphlayouts/reference/layout_focus_group.md)
The layout_igraph\_\* function should not be used directly. It is only
used as an argument for plotting with 'igraph'. 'ggraph' natively
supports the layout.

## Examples

``` r
library(igraph)
g <- sample_gnp(10, 0.4)
coords <- layout_with_focus(g, v = 1)
coords
#> $xy
#>             [,1]        [,2]
#>  [1,] 0.00000000  0.00000000
#>  [2,] 0.92575271 -0.37812951
#>  [3,] 0.96867819  0.24831946
#>  [4,] 1.85779948 -0.74066259
#>  [5,] 1.33574107 -1.48855494
#>  [6,] 2.98499984 -0.29962638
#>  [7,] 0.08349121 -0.99650851
#>  [8,] 0.51811159 -0.85531303
#>  [9,] 1.51507731 -1.30558062
#> [10,] 1.99806892 -0.08786698
#> 
#> $distance
#>  [1] 0 1 1 2 2 3 1 1 2 2
#> 
```
