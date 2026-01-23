# Changelog

## graphlayouts 1.2.2.9000

- use air formatter
- refactoring C++ code regarding stress majorization

## graphlayouts 1.2.2

CRAN release: 2025-01-23

- fixed a bug in multilevel layouts that prevented proper handling of
  layouts in lists

## graphlayouts 1.2.1

CRAN release: 2024-11-18

- moves oaqc back to suggested packages and removed ported code

## graphlayouts 1.2.0

CRAN release: 2024-09-24

- ported relevant code from archived R package oaqc
  ([\#83](https://github.com/schochastics/graphlayouts/issues/83))
- fixed igraph deprecation warnings and require igraph \>= 2.0.0
- removed vignette and point to tutorial
- removed dependency of ggraph

## graphlayouts 1.1.1

CRAN release: 2024-03-09

- fixed bug in disconnected layouts
  [\#80](https://github.com/schochastics/graphlayouts/issues/80)

## graphlayouts 1.1.0

CRAN release: 2024-01-19

- [`layout_with_constrained_stress()`](https://schochastics.github.io/graphlayouts/reference/layout_constrained_stress.md)
  and
  [`layout_with_constrained_stress3D()`](https://schochastics.github.io/graphlayouts/reference/layout_constrained_stress3D.md)
  work for disconnected graphs
- internal code refactoring
- added
  [`layout_as_metromap()`](https://schochastics.github.io/graphlayouts/reference/layout_as_metromap.md)
- added
  [`layout_with_fixed_coords()`](https://schochastics.github.io/graphlayouts/reference/layout_fixed_coords.md)
- removed deprecated igraph calls

## graphlayouts 1.0.2

CRAN release: 2023-11-03

- fixed bug with weighted disconnected graphs
  ([\#71](https://github.com/schochastics/graphlayouts/issues/71)) h/t
  [@gi0na](https://github.com/gi0na)

## graphlayouts 1.0.1

CRAN release: 2023-09-19

- removed the use of `%u%`
  ([\#70](https://github.com/schochastics/graphlayouts/issues/70))

## graphlayouts 1.0.0

CRAN release: 2023-05-01

- added install of oaqc to
  readme([\#52](https://github.com/schochastics/graphlayouts/issues/52))
- fixed grammar in Description
  ([\#53](https://github.com/schochastics/graphlayouts/issues/53))
- made dynamic layout example reproducible
  ([\#54](https://github.com/schochastics/graphlayouts/issues/54))
- replaced 1:length with seq_along
  ([\#55](https://github.com/schochastics/graphlayouts/issues/55))
- added contributing guide
  ([\#56](https://github.com/schochastics/graphlayouts/issues/56))
- added more tests
  ([\#60](https://github.com/schochastics/graphlayouts/issues/60))

## graphlayouts 0.8.4

CRAN release: 2022-11-24

- added more unit tests
- added package level description

## graphlayouts 0.8.3

CRAN release: 2022-10-20

- fixed error for disconnected graphs with an explicit weights vector
  ([\#47](https://github.com/schochastics/graphlayouts/issues/47))
- added proper citation

## graphlayouts 0.8.2

CRAN release: 2022-09-29

- fixed error for very large graphs
  ([\#45](https://github.com/schochastics/graphlayouts/issues/45))
- added
  [`layout_with_focus_group()`](https://schochastics.github.io/graphlayouts/reference/layout_focus_group.md)
  and
  [`layout_with_centrality_group()`](https://schochastics.github.io/graphlayouts/reference/layout_centrality_group.md)
  ([\#46](https://github.com/schochastics/graphlayouts/issues/46))

## graphlayouts 0.8.1

CRAN release: 2022-08-11

- added warning in
  [`layout_as_backbone()`](https://schochastics.github.io/graphlayouts/reference/layout_backbone.md)
  if graph is disconnected and contains isolates

## graphlayouts 0.8.0

CRAN release: 2022-01-03

- added
  [`layout_with_umap()`](https://schochastics.github.io/graphlayouts/reference/layout_umap.md)

## graphlayouts 0.7.2

CRAN release: 2021-11-21

- fixed description of `bbox` in `layout_with_stress`
- fixed bug in `layout_with_stress3D` which only produced a 2D layout

## graphlayouts 0.7.1

CRAN release: 2020-10-26

- restoring old seed after using stress layout

## graphlayouts 0.7.0

CRAN release: 2020-04-25

- added
  [`layout_as_multilevel()`](https://schochastics.github.io/graphlayouts/reference/layout_multilevel.md)
  for multilevel networks
- added
  [`layout_with_stress3D()`](https://schochastics.github.io/graphlayouts/reference/layout_stress3D.md)
  and
  [`layout_with_constrained_stress3D()`](https://schochastics.github.io/graphlayouts/reference/layout_constrained_stress3D.md)
  for 3D layouts
- fixed crash in
  [`layout_as_backbone()`](https://schochastics.github.io/graphlayouts/reference/layout_backbone.md)
  when the graph has loops
  ([\#32](https://github.com/schochastics/graphlayouts/issues/32))

## graphlayouts 0.6.0

CRAN release: 2020-03-09

- added
  [`layout_with_constrained_stress()`](https://schochastics.github.io/graphlayouts/reference/layout_constrained_stress.md)
- added fixed random seed for stress (stress is deterministic and
  produces same layout up to translation/rotation)
- speedup of
  [`layout_with_sparse_stress()`](https://schochastics.github.io/graphlayouts/reference/layout_sparse_stress.md)
  and
  [`layout_with_pmds()`](https://schochastics.github.io/graphlayouts/reference/layout_pmds.md)
  by “smarter” distance calculation
- speedup of
  [`layout_with_sparse_stress()`](https://schochastics.github.io/graphlayouts/reference/layout_sparse_stress.md)
  by using precomputed distances in
  [`layout_with_pmds()`](https://schochastics.github.io/graphlayouts/reference/layout_pmds.md)
- speedup of
  [`layout_with_stress()`](https://schochastics.github.io/graphlayouts/reference/layout_stress.md)
  by dynamically switching to
  [`layout_with_pmds()`](https://schochastics.github.io/graphlayouts/reference/layout_pmds.md)
  during initialisation for large graphs

## graphlayouts 0.5.0

CRAN release: 2019-08-20

- **BREAKING CHANGE**: removed `qgraph()`. Now part of `ggraph`.
- **POSSIBLE BREAKING CHANGE**:
  [`layout_with_focus()`](https://schochastics.github.io/graphlayouts/reference/layout_focus.md)
  now also returns the distance to the focus node
- changed filenames (doesn’t have any effect on functionality)
- added
  [`layout_as_dynamic()`](https://schochastics.github.io/graphlayouts/reference/layout_dynamic.md)
  for longitudinal network data
- removed `gbp`and `scales` dependency and moved `oaqc` to suggest
- edge weights are now supported in
  [`layout_with_stress()`](https://schochastics.github.io/graphlayouts/reference/layout_stress.md)
  and
  [`layout_with_focus()`](https://schochastics.github.io/graphlayouts/reference/layout_focus.md)
- added
  [`layout_with_pmds()`](https://schochastics.github.io/graphlayouts/reference/layout_pmds.md)
  (Pivot MDS for large graphs)
- added
  [`layout_with_sparse_stress()`](https://schochastics.github.io/graphlayouts/reference/layout_sparse_stress.md)
  (“stress for large graphs”)

## graphlayouts 0.2.0

CRAN release: 2019-07-04

- added checks for multiple and directed edges in
  [`layout_as_backbone()`](https://schochastics.github.io/graphlayouts/reference/layout_backbone.md)
- faster implementation of reweighting for
  [`layout_as_backbone()`](https://schochastics.github.io/graphlayouts/reference/layout_backbone.md)
- minor bug fixes in “stress” calculation

## graphlayouts 0.1.0

CRAN release: 2019-04-05

- added more examples and documentation

## graphlayouts 0.0.5.9000

- changed name from `smglr` to `graphlayouts` (sorry)
- added
  [`layout_with_eigen()`](https://schochastics.github.io/graphlayouts/reference/layout_spectral.md)
- layouts can now be used directly in ggraph,
  e.g. `ggraph(g,layout="stress")+...`
- added documentation, examples and references
