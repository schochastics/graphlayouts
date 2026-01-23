# radial focus group layout

arrange nodes in concentric circles around a focal node according to
their distance from the focus and keep predefined groups in the same
angle range.

## Usage

``` r
layout_with_focus_group(
  g,
  v,
  group,
  shrink = 10,
  weights = NA,
  iter = 500,
  tol = 1e-04
)

layout_igraph_focus_group(
  g,
  v,
  group,
  shrink = 10,
  weights = NA,
  iter = 500,
  tol = 1e-04,
  circular
)
```

## Arguments

- g:

  igraph object

- v:

  id of focal node to be placed in the center

- group:

  vector indicating grouping of nodes

- shrink:

  shrink the reserved angle range for a group to increase the gaps
  between groups

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

matrix of xy coordinates

## Details

Be careful when using weights. In most cases, the inverse of the edge
weights should be used to ensure that the endpoints of an edges with
higher weights are closer together (weights=1/E(g)\$weight).

## See also

[layout_focus](https://schochastics.github.io/graphlayouts/reference/layout_focus.md)
The layout_igraph\_\* function should not be used directly. It is only
used as an argument for plotting with 'igraph'.

## Examples

``` r
library(igraph)
g <- sample_islands(4, 5, 0.8, 2)
grp <- as.character(rep(1:4, each = 5))
layout_with_focus_group(g, v = 1, group = grp, shrink = 10)
#>             [,1]       [,2]
#>  [1,]  0.0000000  0.0000000
#>  [2,]  0.8510684  0.5250548
#>  [3,]  0.7695020  0.6386444
#>  [4,]  0.7115721  0.7026131
#>  [5,]  0.1736482  0.9848078
#>  [6,] -1.0855186  1.6797766
#>  [7,] -1.7355044  0.9939942
#>  [8,] -0.3472964  1.9696155
#>  [9,] -1.9696155  0.3472964
#> [10,] -0.5953650  0.8034554
#> [11,] -2.4079115 -1.7894027
#> [12,] -2.9544233 -0.5209445
#> [13,] -1.5238498 -1.2953308
#> [14,] -0.5209445 -2.9544233
#> [15,] -0.3977501 -1.9600497
#> [16,]  2.7205401 -1.2643819
#> [17,]  1.9696155 -0.3472964
#> [18,]  0.5209445 -2.9544233
#> [19,]  1.8589997 -2.3545955
#> [20,]  1.4615800 -1.3652048
```
