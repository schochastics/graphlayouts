# dynamic graph layout

Create layouts for longitudinal networks.

## Usage

``` r
layout_as_dynamic(gList, weights = NA, alpha = 0.5, iter = 500, tol = 1e-04)
```

## Arguments

- gList:

  list of igraph objects. Each network must contain the same set of
  nodes.

- weights:

  possibly a numeric vector with edge weights. If this is NULL and the
  graph has a weight edge attribute, then the attribute is used. If this
  is NA then no weights are used (even if the graph has a weight
  attribute). By default, weights are ignored. See details for more.

- alpha:

  weighting of reference layout. See details.

- iter:

  number of iterations during stress optimization

- tol:

  stopping criterion for stress optimization

## Value

list of coordinates for each graph

## Details

The reference layout is calculated based on the union of all graphs. The
parameter alpha controls the influence of the reference layout. For
alpha=1, only the reference layout is used and all graphs have the same
layout. For alpha=0, the stress layout of each individual graph is used.
Values in-between interpolate between the two layouts.

Be careful when using weights. In most cases, the inverse of the edge
weights should be used to ensure that the endpoints of an edges with
higher weights are closer together (weights=1/E(g)\$weight).

## References

Brandes, U. and Indlekofer, N. and Mader, M. (2012). Visualization
methods for longitudinal social networks and stochastic actor-oriented
modeling. *Social Networks* 34 (3) 291-308

## Examples

``` r
library(igraph)
g1 <- sample_gnp(20, 0.2)
g2 <- sample_gnp(20, 0.2)
g3 <- sample_gnp(20, 0.2)

xy <- layout_as_dynamic(list(g1, g2, g3))

# layout for first network
xy[[1]]
#>             [,1]        [,2]
#>  [1,] -0.1364685  0.03544803
#>  [2,]  0.9608016 -2.02309586
#>  [3,] -1.1069562  1.71803325
#>  [4,] -1.6248427  0.09862508
#>  [5,] -0.5843827 -0.55697364
#>  [6,]  1.7254988  0.12491723
#>  [7,]  0.2333246 -0.69965990
#>  [8,] -0.7571420  0.67668072
#>  [9,]  0.7861468 -0.03389099
#> [10,]  0.1538646 -1.73636243
#> [11,] -0.3395158  1.33874942
#> [12,]  0.4047310  1.42555448
#> [13,] -0.8481309 -1.27250868
#> [14,]  0.1104266 -1.33079786
#> [15,]  1.4410278  1.35643706
#> [16,]  0.7943408  0.33533693
#> [17,]  1.2805972 -0.45651907
#> [18,] -0.5420913  2.95906253
#> [19,] -0.1383375 -0.59079071
#> [20,] -2.2144566 -0.91262040
```
