# graphlayouts: layout algorithms for network visualizations

The package implements several new layout algorithms to visualize
networks. Most are based on the concept of stress majorization. Some
more specific algorithms allow to emphasize hidden group structures in
networks or focus on specific nodes. The package is best used in
conjunction with ggraph.

Some features of the package are:

- [`layout_with_stress()`](https://schochastics.github.io/graphlayouts/reference/layout_stress.md)
  is a state of the art deterministic layout algorithms.

- [`layout_as_backbone()`](https://schochastics.github.io/graphlayouts/reference/layout_backbone.md)
  uncovers hidden group structures (if they exist) by emphasizing
  strongly embedded edges.

- [`layout_with_focus()`](https://schochastics.github.io/graphlayouts/reference/layout_focus.md)
  and
  [`layout_with_centrality()`](https://schochastics.github.io/graphlayouts/reference/layout_centrality.md)
  produce concentric layouts with a focal or most central nodes in the
  center.

- [`layout_with_eigen()`](https://schochastics.github.io/graphlayouts/reference/layout_spectral.md)
  implements some layout algorithms on the basis of eigenvectors

- [`layout_with_sparse_stress()`](https://schochastics.github.io/graphlayouts/reference/layout_sparse_stress.md)
  sparse stress for large graphs

- [`layout_with_pmds()`](https://schochastics.github.io/graphlayouts/reference/layout_pmds.md)
  pivot MDS for large graphs.

- [`layout_as_dynamic()`](https://schochastics.github.io/graphlayouts/reference/layout_dynamic.md)
  for longitudinal network data

A detailed tutorial can be found at
<https://schochastics.github.io/netVizR/>

## See also

Useful links:

- <https://github.com/schochastics/graphlayouts>

- <https://schochastics.github.io/graphlayouts/>

- Report bugs at <https://github.com/schochastics/graphlayouts/issues>

## Author

**Maintainer**: David Schoch <david@schochastics.net>
([ORCID](https://orcid.org/0000-0003-2952-4812))
