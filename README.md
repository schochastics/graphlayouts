<!-- README.md is generated from README.Rmd. Please edit that file -->
smglr
=====

[![Travis build
status](https://travis-ci.org/schochastics/smglr.svg?branch=master)](https://travis-ci.org/schochastics/smglr)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/schochastics/smglr?branch=master&svg=true)](https://ci.appveyor.com/project/schochastics/smglr)

smglr stands for **S**tress **M**ajorization **G**raph **L**ayout in R.
See my dedicated [blog
post](http://blog.schochastics.net/post/stress-based-graph-layouts/) for
more information.

Example: Connected Network
--------------------------

*This example is a bit of a special case and exploits some weird issues
in igraph. Working on adding something more representative*

``` r
library(igraph)   
library(ggraph)   
# devtools::install_github("schochastics/smglr")
library(smglr)

set.seed(666)
pa <- sample_pa(1000,1,1,directed = F)

ggraph(pa)+
  geom_edge_link(width=0.2,colour="grey")+
  geom_node_point(col="black",size=0.3)+
  theme_graph()
```

<img src="figures/README-example-1.png" width="80%" style="display: block; margin: auto;" />

``` r


l <- stress_majorization(pa)
ggraph(pa,layout="manual",node.positions=data.frame(x=l[,1],y=l[,2]))+
  geom_edge_link(width=0.2,colour="grey")+
  geom_node_point(col="black",size=0.3)+
  theme_graph()
```

<img src="figures/README-example-2.png" width="80%" style="display: block; margin: auto;" />

Example: Unconnected Network
----------------------------

Stress majorization now also works for networks with several components
(but not very well yet). It relies on a bin packing algorithm to
efficiently put the components in a rectangle, rather than a circle.

``` r
set.seed(666)
g <- disjoint_union(
  sample_pa(10,directed = F),
  sample_pa(20,directed = F),
  sample_pa(30,directed = F),
  sample_pa(40,directed = F),
  sample_pa(50,directed = F),
  sample_pa(60,directed = F),
  sample_pa(80,directed = F)
)

ggraph(g) +
  geom_edge_link() +
  geom_node_point() +
  theme_graph()
```

<img src="figures/README-example_un-1.png" width="80%" style="display: block; margin: auto;" />

``` r

l <- stress_majorization(g,bbox=30)
ggraph(g, layout="manual", node.positions=data.frame(x=l[,1],y=l[,2])) +
  geom_edge_link() +
  geom_node_point() +
  theme_graph()
```

<img src="figures/README-example_un-2.png" width="80%" style="display: block; margin: auto;" />
