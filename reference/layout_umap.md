# UMAP graph layouts

Using the UMAP dimensionality reduction algorithm as a graph layout

## Usage

``` r
layout_with_umap(g, pivots = NULL, ...)

layout_igraph_umap(g, circular, ...)
```

## Arguments

- g:

  igraph object

- pivots:

  if not NULL, number of pivot nodes to use for distance calculation
  (for large graphs).

- ...:

  additional parameters for umap. See the
  [`?uwot::umap`](https://jlmelville.github.io/uwot/reference/umap.html)
  for help.

- circular:

  not used

## Value

matrix of xy coordinates

## Details

The layout_igraph\_\* function should not be used directly. It is only
used as an argument for plotting with 'igraph'. UMAP can be tuned by
many different parameters. Refer to the documentation at
https://github.com/jlmelville/uwot for help

## References

McInnes, Leland, John Healy, and James Melville. "Umap: Uniform manifold
approximation and projection for dimension reduction." arXiv preprint
arXiv:1802.03426 (2018).

## Author

David Schoch

## Examples

``` r
library(igraph)

g <- sample_islands(10, 20, 0.6, 10)
# xy <- layout_with_umap(g, min_dist = 0.5)
```
