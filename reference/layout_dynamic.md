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
#>  [1,]  0.09352920  1.76120679
#>  [2,]  0.14695919  0.81110145
#>  [3,]  1.87090391  0.51814189
#>  [4,] -1.39186189  0.33500741
#>  [5,]  0.32116974 -1.18268736
#>  [6,] -1.23872085 -0.38963547
#>  [7,] -0.24433495 -0.67069030
#>  [8,]  0.61365422 -0.01501541
#>  [9,] -0.41810078  1.07421005
#> [10,]  0.82747280 -0.10684176
#> [11,]  1.00269008 -1.30356798
#> [12,] -1.27276186  0.19156259
#> [13,]  2.13515963 -0.69744071
#> [14,] -0.07851292 -1.21895194
#> [15,] -1.00289180  1.63518840
#> [16,] -0.19870658 -2.45551354
#> [17,] -0.90463472 -0.55902907
#> [18,] -1.54491881  0.15983749
#> [19,]  0.15647408 -0.47947715
#> [20,]  0.93483943  1.82185256
```
