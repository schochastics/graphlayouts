# multilevel layout

Layout algorithm to visualize multilevel networks

## Usage

``` r
layout_as_multilevel(
  g,
  type = "all",
  FUN1,
  FUN2,
  params1 = NULL,
  params2 = NULL,
  ignore_iso = TRUE,
  project2D = TRUE,
  alpha = 35,
  beta = 45
)

layout_igraph_multilevel(
  g,
  type = "all",
  FUN1,
  FUN2,
  params1 = NULL,
  params2 = NULL,
  ignore_iso = TRUE,
  alpha = 35,
  beta = 45,
  circular
)
```

## Arguments

- g:

  igraph object. Must have a vertex attribute "lvl" which is 1 or 2.

- type:

  one of "all", "separate","fix1" or "fix2". see details

- FUN1:

  if type="separate", the layout function to be used for level 1

- FUN2:

  if type="separate", the layout function to be used for level 2

- params1:

  named list of parameters for FUN1

- params2:

  named list of parameters for FUN2

- ignore_iso:

  treatment of isolates within levels. see details

- project2D:

  logical. Defaults to TRUE (project to 2D).

- alpha:

  angle for isometric projection between 0 and 90

- beta:

  angle for isometric projection between 0 and 90

- circular:

  not used

## Value

matrix of xy coordinates

## Details

The algorithm internally computes a 3D layout where each level is in a
separate y-plane. The layout is then projected into 2D via an isometric
mapping, controlled by the parameters `alpha` and `beta`. It may take
some adjusting to `alpha` and `beta` to find a good perspective.

If type="all", the layout is computed at once for the complete network.
For type="separate", two user specified layout algorithms (`FUN1` and
`FUN2`) are used for the levels. The named lists `param1` and `param2`
can be used to set parameters for `FUN1` and `FUN2`. This option helpful
for situations where different structural features of the levels should
be emphasized.

For type="fix1" and type="fix2" only one of the level layouts is fixed.
The other one is calculated by optimizing the inter level ties, such
that they are drawn (almost) vertical.

The `ignore_iso` parameter controls the handling of isolates. If TRUE,
nodes without inter level edges are ignored during the layout process
and added at the end. If FALSE they are left unchanged

The layout_igraph\_\* function should not be used directly. It is only
used as an argument for plotting with 'igraph'.

## Examples

``` r
library(igraph)
data("multilvl_ex")
if (FALSE) { # \dontrun{
# compute a layout for the whole network
xy <- layout_as_multilevel(multilvl_ex, type = "all", alpha = 25, beta = 45)

# compute a layout for each level separately and combine them
xy <- layout_as_multilevel(multilvl_ex,
    type = "separate",
    FUN1 = layout_as_backbone,
    FUN2 = layout_with_stress,
    alpha = 25, beta = 45
)
} # }
```
