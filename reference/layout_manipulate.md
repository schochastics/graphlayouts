# manipulate layout

functions to manipulate an existing layout

## Usage

``` r
layout_rotate(xy, angle)

layout_mirror(xy, axis = "vertical")
```

## Arguments

- xy:

  graph layout

- angle:

  angle for rotation

- axis:

  mirror horizontal or vertical

## Value

manipulated matrix of xy coordinates

## Details

These functions are mostly useful for deterministic layouts such as
[layout_with_stress](https://schochastics.github.io/graphlayouts/reference/layout_stress.md)

## Author

David Schoch

## Examples

``` r
library(igraph)
g <- sample_gnp(50, 0.3)

xy <- layout_with_stress(g)

# rotate 90 degrees
xy <- layout_rotate(xy, 90)

# flip horizontally
xy <- layout_mirror(xy, "horizontal")
```
