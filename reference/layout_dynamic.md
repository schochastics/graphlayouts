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
#>              [,1]        [,2]
#>  [1,]  0.09331900  1.76148108
#>  [2,]  0.14873289  0.81045551
#>  [3,]  1.87091946  0.51753695
#>  [4,] -1.39373120  0.33579754
#>  [5,]  0.32103219 -1.18210059
#>  [6,] -1.23820495 -0.38965935
#>  [7,] -0.24221774 -0.67195390
#>  [8,]  0.61201479 -0.01115133
#>  [9,] -0.41951089  1.07427409
#> [10,]  0.82694504 -0.10531909
#> [11,]  1.00149105 -1.30723378
#> [12,] -1.27220013  0.19152664
#> [13,]  2.13481069 -0.69603347
#> [14,] -0.07613157 -1.21936832
#> [15,] -1.00379044  1.63416677
#> [16,] -0.19754008 -2.45558780
#> [17,] -0.90351839 -0.55927490
#> [18,] -1.54436579  0.15647447
#> [19,]  0.15502943 -0.47713634
#> [20,]  0.93345936  1.82275581
```
