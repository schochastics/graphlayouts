---
title: 'graphlayouts: Layout algorithms for network visualizations in R'
tags:
- R
- network visualization
- graph drawing
- layout algorithms
date: "03 March 2022"
output: pdf_document
affiliations:
- name: "GESIS - Leibniz Institute for the Social Sciences"
  index: 1
authors:
- name: David Schoch
  orcid: "0000-0003-2952-4812"
  affiliation: 1
bibliography: paper.bib
---


# Summary
Network analysis is a powerful methodology to understand the relationships among entities, such as people, or organizations, and the patterns that emerge from these connections. It is increasingly popular in various fields, including sociology, biology, economics, and computer science, and has been used to study diverse phenomena including the spread of diseases, flow of information, and the structure of political organizations. Network visualization is a powerful tool for exploring, analyzing, and communicating network structures and patterns therein. However, in contrast to tabular data, nodes can technically be placed arbitrarily on the plane and it is easy to draw wrong conclusion based on an inadequate layout. To circumvent arbitrary placement of nodes, many different layout algorithms have been developed which optimize different stylistic features and can serve purposes in communicating structural properties. The package `graphlayouts` implements several state-of-the-art algorithms which are so far not available for R. This includes algorithms for large graphs, to emphasize hidden group structures, and important nodes within a network.

# Statement of need

The all-purpose network analysis package `igraph` [@cn-ispcnr-06] already implements a 
great variety of layout algorithms for networks. `graphlayouts` complements these by adding a faster general purpose
algorithm and a series of specialized algorithms that serve very specific purposes, either to emphasize group structures or
the position of individual nodes within a network. It further adds support for dynamic and multilevel networks.
The package is already well integrated into the R ecosystem. `ggraph` [@p-g-22], the `ggplot2` for networks, imports the package and uses the stress based layout as its default layout algorithm. All layout algorithms included in the package by default return a matrix of coordinates. This allows to use the layouts with most other network visualization packages too, including static ones such as `sna` [@b-snas-08] and `statnet` [@handcock2008statnet], and interactive ones such as `visNetwork` [@visNetwork2022].

# Overview of algorithms

In this section, the most prominent implemented layout algorithms are introduced. Others
can be found in the online documentation linked above.

More details on plotting with `graphlayouts` can be found in a short vignette, the online
documentation <http://graphlayouts.schochastics.net> and a dedicated tutorial on
network visualization <https://www.mr.schochastics.net/material/netVizR/>.

## Stress majorization

Stress majorization [@gansner2005graph] is, in contrast to many other layout
algorithms, deterministic and able to produce high quality layouts for a great
variety of graph classes. `graphlayouts` implements classic stress majorization
in `layout_with_stress()`. An example is shown in Figure \ref{fig:stress-ex}.

```r
ggraph(pa,layout = "stress")+
  geom_edge_link0(edge_width = 0.2, edge_colour = "grey")+
  geom_node_point(size = 0.3)+
  theme_graph()
```

\begin{figure}[htb]
  \centering
  \includegraphics[width=0.75\textwidth]{man/figures/README-example-2.png}
  \caption{Example of a stress based layout.}
  \label{fig:stress-ex}
\end{figure}

## Sparse stress majorization

Stress majorization requires the computation of the full distance matrix, which becomes
computationally expensive for large graphs. The function
`layout_with_sparse_stress()` calculates the distances only for a small set of
pivots and creates the layout based on these distances[@ortmann2016sparse]. The resulting layout is
not as good as the full stress, but it can reasonably be used for graphs with
around 100k nodes and 5 million edges. The wiki
<https://github.com/schochastics/graphlayouts/wiki> contains several benchmark
results in comparison with layouts from `igraph`.

## Backbone layout

`layout_as_backbone()` is a layout algorithm that can help emphasize hidden
group structures in hairball graphs [@nocaj2015untangling]. To illustrate the performance of the
algorithm, we use an
artificial network with a subtle group structure (cf. Figure \ref{fig:hairball}).
```r
set.seed(665)
#create network with a group structure
g <- sample_islands(9,40,0.4,15)
g <- simplify(g)
V(g)$grp <- as.character(rep(1:9,each=40))

ggraph(g,layout = "stress")+
  geom_edge_link0(colour = rgb(0,0,0,0.5), edge_width = 0.1)+
  geom_node_point(aes(col = grp))+
  scale_color_brewer(palette = "Set1")+
  theme_graph()+
  theme(legend.position = "none")
```
\begin{figure}
  \centering
  \includegraphics[width=0.75\textwidth]{man/figures/README-hairball-1.png}
  \caption{Hairball graph with stress layout.}
  \label{fig:hairball}
\end{figure}

The backbone layout helps to uncover potential group structures based on edge
embeddedness and puts more emphasis on this structure in the layout (cf. Figure \ref{fig:backbone}).

```r
bb <- layout_as_backbone(g, keep = 0.4)
E(g)$col <- FALSE
E(g)$col[bb$backbone] <- TRUE

ggraph(g, layout = "manual", x = bb$xy[, 1], y = bb$xy[, 2]) +
  geom_edge_link0(aes(col = col), edge_width = 0.1) +
  geom_node_point(aes(col = grp)) +
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
  theme_graph() +
  theme(legend.position = "none")
```

\begin{figure}
  \centering
  \includegraphics[width=0.75\textwidth]{man/figures/README-backbone-1.png}
  \caption{Backbone layout of the graph shown in Figure \ref{fig:hairball}.}
  \label{fig:backbone}
\end{figure}

## Radial layouts

The function `layout_with_focus()` creates a radial layout around a
focal node [@brandes2010more]. All nodes with the same distance from the focal node are on
the same circle (cf Figure \ref{fig:focus}).

``` r
library(igraphdata)
library(patchwork)
data("karate")

p1 <- ggraph(karate, layout = "focus", focus = 1) +
  draw_circle(use = "focus", max.circle = 3) +
  geom_edge_link0(edge_color = "black", edge_width = 0.3) +
  geom_node_point(aes(fill = as.factor(Faction)), size = 2, shape = 21) +
  scale_fill_manual(values = c("#8B2323", "#EEAD0E")) +
  theme_graph() +
  theme(legend.position = "none") +
  coord_fixed() +
  labs(title = "Focus on Mr. Hi")

p2 <- ggraph(karate, layout = "focus", focus = 34) +
  draw_circle(use = "focus", max.circle = 4) +
  geom_edge_link0(edge_color = "black", edge_width = 0.3) +
  geom_node_point(aes(fill = as.factor(Faction)), size = 2, shape = 21) +
  scale_fill_manual(values = c("#8B2323", "#EEAD0E")) +
  theme_graph() +
  theme(legend.position = "none") +
  coord_fixed() +
  labs(title = "Focus on John A.")

p1 + p2
```

\begin{figure}
  \centering
  \includegraphics[width=0.8\textwidth]{man/figures/README-flex_focus-1.png}
  \caption{Example of focus layouts.}
  \label{fig:focus}
\end{figure}

The function `layout_with_centrality` creates a radial layout around the
node with the highest centrality value. The further outside a node is,
the more peripheral it is (cf. Figure \ref{fig:cent}).

``` r
bc <- betweenness(karate)
p1 <- ggraph(karate, layout = "centrality", centrality = bc, tseq = seq(0, 1, 0.15)) +
  draw_circle(use = "cent") +
  annotate_circle(bc, format = "", pos = "bottom") +
  geom_edge_link0(edge_color = "black", edge_width = 0.3) +
  geom_node_point(aes(fill = as.factor(Faction)), size = 2, shape = 21) +
  scale_fill_manual(values = c("#8B2323", "#EEAD0E")) +
  theme_graph() +
  theme(legend.position = "none") +
  coord_fixed() +
  labs(title = "betweenness centrality")


cc <- closeness(karate)
p2 <- ggraph(karate, layout = "centrality", centrality = cc, tseq = seq(0, 1, 0.2)) +
  draw_circle(use = "cent") +
  annotate_circle(cc, format = "scientific", pos = "bottom") +
  geom_edge_link0(edge_color = "black", edge_width = 0.3) +
  geom_node_point(aes(fill = as.factor(Faction)), size = 2, shape = 21) +
  scale_fill_manual(values = c("#8B2323", "#EEAD0E")) +
  theme_graph() +
  theme(legend.position = "none") +
  coord_fixed() +
  labs(title = "closeness centrality")

p1 + p2
```

\begin{figure}
  \centering
  \includegraphics[width=0.8\textwidth]{man/figures/README-flex_cent-1.png}
  \caption{Example of centrality layout using betweenness and closeness.}
  \label{fig:cent}
\end{figure}

# References
