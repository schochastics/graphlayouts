test_that("layout_dynamic works", {
  library(igraph)
  expect_error(layout_as_dynamic(3))
  expect_error(layout_as_dynamic(list(igraph::graph.full(3),igraph::graph.full(5))))
  expect_is(layout_as_dynamic(list(igraph::graph.full(5),igraph::graph.full(5))),"list")
})
