# Manipulate graph

functions to manipulate a graph

## Usage

``` r
reorder_edges(g, attr, desc = TRUE)
```

## Arguments

- g:

  igraph object

- attr:

  edge attribute name used to sort edges

- desc:

  logical. sort in descending (default) or ascending order

## Value

manipulated graph

## Details

`reorder_edges()` allows to reorder edges according to an attribute so
that edges are drawn in the given order.

## Author

David Schoch

## Examples

``` r
library(igraph)

g <- sample_gnp(10, 0.5)
E(g)$attr <- 1:ecount(g)
gn <- reorder_edges(g,"attr")
```
