# annotate concentric circles

annotate concentric circles

## Usage

``` r
annotate_circle(cent, col = "#00BFFF", format = "", pos = "top", text_size = 3)
```

## Arguments

- cent:

  centrality scores used for layout

- col:

  color of text

- format:

  either empty string or 'scientific'

- pos:

  position of text ('top' or 'bottom')

- text_size:

  font size for annotations

## Value

annotated concentric circles around origin

## Details

this function is best used with
[layout_with_centrality](https://schochastics.github.io/graphlayouts/reference/layout_centrality.md)
together with
[draw_circle](https://schochastics.github.io/graphlayouts/reference/draw_circle.md).

## Examples

``` r
library(igraph)
#> 
#> Attaching package: ‘igraph’
#> The following objects are masked from ‘package:stats’:
#> 
#>     decompose, spectrum
#> The following object is masked from ‘package:base’:
#> 
#>     union

g <- sample_gnp(10, 0.4)
if (FALSE) { # \dontrun{
library(ggraph)
ggraph(g, layout = "centrality", centrality = closeness(g)) +
    draw_circle(use = "cent") +
    annotate_circle(closeness(g), pos = "bottom", format = "scientific") +
    geom_edge_link() +
    geom_node_point(shape = 21, fill = "grey25", size = 5) +
    theme_graph() +
    coord_fixed()
} # }
```
