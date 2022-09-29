# graphlayouts 0.8.2

* fixed error for very large graphs (#45)
* added `layout_with_focus_group()` and `layout_with_centrality_group()` (#46)

# graphlayouts 0.8.1

* added warning in `layout_as_backbone()` if graph is disconnected and contains isolates

# graphlayouts 0.8.0

* added `layout_with_umap()`

# graphlayouts 0.7.2

* fixed description of `bbox` in `layout_with_stress`
* fixed bug in `layout_with_stress3D` which only produced a 2D layout

# graphlayouts 0.7.1

* restoring old seed after using stress layout

# graphlayouts 0.7.0

* added `layout_as_multilevel()` for multilevel networks
* added `layout_with_stress3D()` and `layout_with_constrained_stress3D()` for 3D layouts
* fixed  crash in `layout_as_backbone()` when the graph has loops (#32)

# graphlayouts 0.6.0

* added `layout_with_constrained_stress()`
* added fixed random seed for stress (stress is deterministic and produces same layout up to translation/rotation)
* speedup of `layout_with_sparse_stress()` and `layout_with_pmds()` by "smarter" distance calculation
* speedup of `layout_with_sparse_stress()` by using precomputed distances in `layout_with_pmds()`
* speedup of `layout_with_stress()` by dynamically switching to `layout_with_pmds()` during initialisation for large graphs


# graphlayouts 0.5.0

* **BREAKING CHANGE**: removed `qgraph()`. Now part of `ggraph`.
* **POSSIBLE BREAKING CHANGE**: `layout_with_focus()` now also returns the distance to the focus node
* changed filenames (doesn't have any effect on functionality)
* added `layout_as_dynamic()` for longitudinal network data
* removed `gbp`and `scales` dependency and moved `oaqc` to suggest
* edge weights are now supported in `layout_with_stress()` and `layout_with_focus()`
* added `layout_with_pmds()` (Pivot MDS for large graphs)
* added `layout_with_sparse_stress()` ("stress for large graphs")

# graphlayouts 0.2.0

* added checks for multiple and directed edges in `layout_as_backbone()`
* faster implementation of reweighting for `layout_as_backbone()`
* minor bug fixes in "stress" calculation

# graphlayouts 0.1.0

* added more examples and documentation

# graphlayouts 0.0.5.9000

* changed name from `smglr` to `graphlayouts` (sorry)
* added `layout_with_eigen()`
* layouts can now be used directly in ggraph, e.g. `ggraph(g,layout="stress")+...`
* added documentation, examples and references

# smglr 0.0.4.9000

* added `layout_with_centrality()`
* added `layout_with_focus()`
* added `reorder_edges()` to reorder an edgelist to control the order of plotting edges

# smglr 0.0.3.9000

* added `layout_as_backbone()`

# smglr 0.0.2.9000

* added functions to rotate/mirror layouts
* `stress_majorization()` is now deprecated. Use `layout_with_stress()` instead

# smglr 0.0.1.9000

* added support for unconnected networks via external library

# smglr 0.0.0.9000

* added a `NEWS.md` file to track changes to the package.
* basic implementation


