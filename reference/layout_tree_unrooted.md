# unrooted tree layout

arranges the nodes of a tree (or forest) in the plane without a
designated root, using algorithms developed for drawing phylogenies.

## Usage

``` r
layout_as_tree_unrooted(
  g,
  mode = c("equalangle", "equaldaylight", "stress"),
  weights = NA,
  daylight_iter = 15,
  iter = 500,
  tol = 1e-04,
  bbox = 30
)

layout_igraph_tree_unrooted(
  g,
  mode = c("equalangle", "equaldaylight", "stress"),
  weights = NA,
  daylight_iter = 15,
  iter = 500,
  tol = 1e-04,
  bbox = 30,
  circular
)
```

## Arguments

- g:

  igraph object. Each connected component must be a tree (acyclic).

- mode:

  which algorithm to use. One of `"equalangle"` (default),
  `"equaldaylight"` or `"stress"`. See details.

- weights:

  possibly a numeric vector with edge weights, interpreted as branch
  lengths. If this is NULL and the graph has a weight edge attribute,
  then the attribute is used. If this is NA then unit branch lengths are
  used (even if the graph has a weight attribute). By default, unit
  branch lengths are used.

- daylight_iter:

  number of refinement sweeps for `mode = "equaldaylight"`.

- iter:

  number of iterations during stress optimization (only used for
  `mode = "stress"`).

- tol:

  stopping criterion for stress optimization (only used for
  `mode = "stress"`).

- bbox:

  width of layout. Only relevant to determine the placement of
  disconnected components (forests).

- circular:

  not used

## Value

matrix of xy coordinates

## Details

Three algorithms are available:

- `"equalangle"` (Felsenstein, 1989): each subtree is allotted an
  angular wedge proportional to its number of leaves. Linear time, but
  can look crowded on unbalanced trees.

- `"equaldaylight"` (Felsenstein): iteratively equalizes the angular
  gaps ("daylight") around each internal node, starting from the
  equal-angle layout, for a more balanced drawing.

- `"stress"`: stress majorization using the patristic (tree path-length)
  distances as target distances, reusing the package's stress engine.

Branch lengths (edge weights) are respected when present; otherwise unit
lengths are used.

The layout_igraph\_\* function should not be used directly. It is only
used as an argument for plotting with 'igraph'. 'ggraph' natively
supports the layout.

## References

Felsenstein, J. (1989). PHYLIP - Phylogeny Inference Package (Version
3.2). *Cladistics*, 5, 164-166.

Bachmaier, C., Brandes, U., & Schlieper, B. (2005). Drawing phylogenetic
trees. *In International Symposium on Algorithms and Computation* (pp.
1110-1121). Springer, Berlin, Heidelberg.

## Examples

``` r
library(igraph)
g <- make_tree(20, 3, mode = "undirected")

xy <- layout_as_tree_unrooted(g)
xy <- layout_as_tree_unrooted(g, mode = "equaldaylight")
```
