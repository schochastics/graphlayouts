# radial centrality group layout

arranges nodes in concentric circles according to a centrality index and
keeping groups within a angle range

## Usage

``` r
layout_with_centrality_group(g, cent, group, shrink = 10, ...)

layout_igraph_centrality_group(g, cent, group, shrink = 10, circular, ...)
```

## Arguments

- g:

  igraph object

- cent:

  centrality scores

- group:

  vector indicating grouping of nodes

- shrink:

  shrink the reserved angle range for a group to increase the gaps
  between groups

- ...:

  additional arguments to layout_with_centrality The layout_igraph\_\*
  function should not be used directly. It is only used as an argument
  for plotting with 'igraph'. 'ggraph' natively supports the layout.

- circular:

  not used

## Value

matrix of xy coordinates

## See also

[layout_centrality](https://schochastics.github.io/graphlayouts/reference/layout_centrality.md)

## Examples

``` r
library(igraph)
```
