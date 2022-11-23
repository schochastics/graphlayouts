library(igraph)
test_that("layout_spectral works", {
  g <- igraph::make_graph( ~ a -- b -- c -- d:e:f)

  expect_is(layout_with_eigen(g,type = "adjacency",ev = "largest"),"matrix")
  expect_is(layout_with_eigen(g,type = "adjacency",ev = "smallest"),"matrix")
  expect_is(layout_with_eigen(g,type = "laplacian",ev = "largest"),"matrix")
  expect_is(layout_with_eigen(g,type = "laplacian",ev = "smallest"),"matrix")

  expect_warning(layout_with_eigen(igraph::as.directed(g)))
  expect_error(layout_with_eigen(igraph::graph.empty(10,directed = FALSE)))
})
