# spectral graph layouts

Using a set of eigenvectors of matrices associated with a graph as
coordinates

## Usage

``` r
layout_with_eigen(g, type = "laplacian", ev = "smallest")

layout_igraph_eigen(g, type = "laplacian", ev = "smallest", circular)
```

## Arguments

- g:

  igraph object

- type:

  matrix to be used for spectral decomposition. either 'adjacency' or
  'laplacian'

- ev:

  eigenvectors to be used. Either 'smallest' or 'largest'.

- circular:

  not used

## Value

matrix of xy coordinates

## Details

The layout_igraph\_\* function should not be used directly. It is only
used as an argument for plotting with 'igraph'. 'ggraph' natively
supports the layout.

## Author

David Schoch

## Examples

``` r
library(igraph)

g <- sample_gnp(50, 0.2)

xy <- layout_with_eigen(g, type = "adjacency", ev = "largest")

xy <- layout_with_eigen(g, type = "adjacency", ev = "smallest")

xy <- layout_with_eigen(g, type = "laplacian", ev = "largest")

xy <- layout_with_eigen(g, type = "laplacian", ev = "smallest")
```
