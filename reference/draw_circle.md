# Draw concentric circles

Draw concentric circles

## Usage

``` r
draw_circle(col = "#00BFFF", use = "focus", max.circle)
```

## Arguments

- col:

  color of circles

- use:

  one of 'focus' or 'cent'

- max.circle:

  if use = 'focus' specifies the number of circles to draw

## Value

concentric circles around origin

## Details

this function is best used with a concentric layout such as
[layout_with_focus](https://schochastics.github.io/graphlayouts/reference/layout_focus.md)
and
[layout_with_centrality](https://schochastics.github.io/graphlayouts/reference/layout_centrality.md).

## Examples

``` r
library(igraph)
g <- sample_gnp(10, 0.4)

if (FALSE) { # \dontrun{
library(ggraph)
ggraph(g, layout = "centrality", centrality = degree(g)) +
    draw_circle(use = "cent") +
    geom_edge_link() +
    geom_node_point(shape = 21, fill = "grey25", size = 5) +
    theme_graph() +
    coord_fixed()
} # }
```
